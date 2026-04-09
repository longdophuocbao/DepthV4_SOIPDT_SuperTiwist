/**
 * ====================================================================
 *  Tractor Depth Controller – ESP32
 * ====================================================================
 *  Mô hình đối tượng (SOIPDT):
 *    G(s) = K · e^(-Ls) / (s(τ₁s + 1))
 *    Phương trình trạng thái: α̈ = −a1·α̇ + b·u
 *      a1 = 1/τ₁,  a2 = 0,  b = K/τ₁
 *
 *  Bù dead-time: Smith Predictor
 *    y_pred(t) = y_measured(t) + y_m(t) − y_m_d(t)
 *    (y_m: model không trễ, y_m_d: model có trễ L)
 *
 *  Mặt trượt PID:
 *    s = Kp·e + Ki·∫e·dt + Kd·ė
 *
 *  Luật điều khiển tổng: u = u_eq + u_sw
 *
 *  Lực tương đương (u_eq) – giải ṡ = 0:
 *    u_eq = (1/b)·(a1·α̇_ctrl + a2·α_ctrl + α̈_ref − (Kp/Kd)·ė − (Ki/Kd)·e)
 *
 *  Lực chuyển mạch (u_sw) – Super-Twisting Algorithm:
 *    u_sw = (1/(b·Kd))·(−K1·√|s|·sgn(s) + ∫−K2·sgn(s)·dt)
 *
 *  Board    : ESP32 Dev Module  (PlatformIO: esp32doit-devkit-v1)
 * ====================================================================
 */

#define DEBUG_SERIAL // tuning by serial
#define MONITOR_WIFI // monitor by wifi

#include <Arduino.h>

#ifdef MONITOR_WIFI
#include <WiFi.h>
#include <LittleFS.h>
#include <AsyncTCP.h>
#include <ESPAsyncWebServer.h>
#include <DNSServer.h>
#endif

#include <math.h>

#include <Wire.h>
#include <Adafruit_ADS1X15.h>
#include <ESP32Servo.h>

// I2C cho ADS1115
#define PIN_SDA 25
#define PIN_SCL 26

// H-bridge: 2 kênh PWM độc lập
#define PIN_MOTOR_IN1 18 // chân PWM kênh thuận
#define PIN_MOTOR_IN2 19 // chân PWM kênh nghịch

#define PWM_FREQ_HZ 20000
#define PWM_BITS 10
#define PWM_MAX_ABS 0.85f

Adafruit_ADS1115 liftingsensor;
Adafruit_ADS1115 tailboardsensor;
ESP32PWM pwmIN1;
ESP32PWM pwmIN2;

#ifdef MONITOR_WIFI
static const char *AP_SSID = "TractorControl";
static const char *AP_PASS = "12345678";
#endif

// ────────────────────────────────────────────────────────────────────
//  CHU KỲ LẤY MẪU
// ────────────────────────────────────────────────────────────────────
static const float DT = 0.002f; // 2 ms = 500 Hz
static const TickType_t DT_TICKS = pdMS_TO_TICKS(2);

static const TickType_t SENSOR_DT_TICKS = pdMS_TO_TICKS(2);

#ifdef MONITOR_WIFI
static const TickType_t WEB_DT_TICKS = pdMS_TO_TICKS(50);
#endif

#ifdef DEBUG_SERIAL
static const TickType_t SERIALDEBUG_DT_TICKS = pdMS_TO_TICKS(20); // 20ms =50hz
#endif

// ────────────────────────────────────────────────────────────────────
//  CẤU TRÚC MÔ HÌNH SOIPDT (Discrete-time, Euler forward)
// ────────────────────────────────────────────────────────────────────
static const int MAX_DELAY_SAMPLES = 2000;

struct SOIPDTModel
{
    float K = 1.0f;
    float tau1 = 5.0f;
    float L = 1.0f;

    float y_k1 = 0.0f;
    float y_k2 = 0.0f;

    float u_buf[MAX_DELAY_SAMPLES];
    int buf_head = 0;
    int d = 0;

    void init(float dt)
    {
        d = (int)roundf(L / dt);
        if (d >= MAX_DELAY_SAMPLES)
            d = MAX_DELAY_SAMPLES - 1;
        memset(u_buf, 0, sizeof(u_buf));
        y_k1 = y_k2 = 0.0f;
        buf_head = 0;
    }

    float step(float u_in, float dt)
    {
        u_buf[buf_head] = u_in;
        int idx = (buf_head - d + MAX_DELAY_SAMPLES) % MAX_DELAY_SAMPLES;
        float u_d = u_buf[idx];
        buf_head = (buf_head + 1) % MAX_DELAY_SAMPLES;

        float dt2 = dt * dt;
        float denom = tau1 + dt;
        float y_n = (K * dt2 * u_d + (2.0f * tau1 + dt) * y_k1 - tau1 * y_k2) / denom;

        y_k2 = y_k1;
        y_k1 = y_n;
        return y_n;
    }

    void reset()
    {
        y_k1 = y_k2 = 0.0f;
        buf_head = 0;
        memset(u_buf, 0, sizeof(u_buf));
    }
};

// ────────────────────────────────────────────────────────────────────
//  REFERENCE FILTER (BỘ LỌC TÍN HIỆU ĐẶT BẬC 2)
// ────────────────────────────────────────────────────────────────────
struct ReferenceFilter
{
    float x1, x2;
    float omega;
    float zeta;

    void init(float _omega, float _zeta)
    {
        omega = _omega;
        zeta = _zeta;
        x1 = x2 = 0.0f;
    }

    struct FilterOutput
    {
        float val, dot, ddot;
    };
    FilterOutput update(float r, float dt)
    {
        float error = r - x1;
        float x_ddot = (omega * omega) * error - (2.0f * zeta * omega) * x2;
        x2 += x_ddot * dt;
        x1 += x2 * dt;
        return {x1, x2, x_ddot};
    }
};

ReferenceFilter refFilter;

// ────────────────────────────────────────────────────────────────────
//  BỘ LỌC THÔNG THẤP BẬC 2 (2nd Order LPF - Butterworth)
// ────────────────────────────────────────────────────────────────────
struct LPF2ndOrder
{
    float b0, b1, b2, a1, a2;
    float x1, x2, y1, y2;

    void init(float fc, float fs)
    {
        float w0 = 2.0f * M_PI * fc / fs;
        float cosw0 = cosf(w0);
        float alpha = sinf(w0) / 1.4142f;
        float a0 = 1.0f + alpha;
        b0 = ((1.0f - cosw0) / 2.0f) / a0;
        b1 = (1.0f - cosw0) / a0;
        b2 = ((1.0f - cosw0) / 2.0f) / a0;
        a1 = (-2.0f * cosw0) / a0;
        a2 = (1.0f - alpha) / a0;
        x1 = x2 = y1 = y2 = 0.0f;
    }

    float update(float x)
    {
        float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x;
        y2 = y1;
        y1 = y;
        return y;
    }
};

LPF2ndOrder filterLifting;
LPF2ndOrder filterTailboard;

// ────────────────────────────────────────────────────────────────────
//  BIẾN TOÀN CỤC – TRẠNG THÁI HỆ THỐNG
// ────────────────────────────────────────────────────────────────────
volatile float g_liftingangle = 40.0f;
volatile float g_tailboardangle = 0.0f;
volatile float g_lifting_filtered = 40.0f;
volatile float g_tail_filtered = 0.0f;
volatile float g_depth_target = 100.0f;
volatile float g_setpoint = 20.0f;
volatile bool g_is_auto = true;
volatile float g_alpha_manual = 20.0f;
volatile float g_e = 0.0f;
volatile float g_de = 0.0f;
volatile float g_s = 0.0f;
volatile float g_u = 0.0f;
volatile float g_u_eq = 0.0f;
volatile float g_u_sw = 0.0f;
volatile float g_u_sw_i = 0.0f;
volatile float g_estimate_depth = 0.0f;
volatile float g_dot_alpha_ref = 0.0f;
volatile float g_ddot_alpha_ref = 0.0f;

// ────────────────────────────────────────────────────────────────────
//  THAM SỐ ĐIỀU KHIỂN
// ────────────────────────────────────────────────────────────────────
volatile float g_Kp = 2.8f;
volatile float g_Ki = 4.0f;
volatile float g_Kd = 1.0f;
volatile float g_K1 = 5.0f;  // Gain của Super-Twisting SMC (√|s|·sign(s))
volatile float g_K2 = 10.0f; // Gain của Super-Twisting SMC (∫sign(s))

volatile float g_K = 1.0f;
volatile float g_tau1 = 5.0f;
volatile float g_L = 1.0f;
volatile float g_lift_offset = 150.0f;
volatile float g_tail_offset = 150.0f;

volatile bool g_model_dirty = true;

// ────────────────────────────────────────────────────────────────────
//  ĐỐI TƯỢNG TOÀN CỤC
// ────────────────────────────────────────────────────────────────────
SemaphoreHandle_t g_mutex;
SOIPDTModel g_model_nd;
SOIPDTModel g_model_wd;

#ifdef MONITOR_WIFI
AsyncWebServer server(80);
AsyncWebSocket ws("/ws");
DNSServer dnsServer;
#endif

float g_lifting_raw_val = 0, g_tail_raw_val = 0, g_sp_raw_val = 0;

uint8_t calculateChecksum(uint8_t *data, size_t len)
{
    uint8_t checksum = 0;
    for (size_t i = 0; i < len; i++)
        checksum ^= data[i];
    return checksum;
}

// ────────────────────────────────────────────────────────────────────
//  HÀM ĐỌC CẢM BIẾN
// ────────────────────────────────────────────────────────────────────
float readLiftingSensorRaw()
{
    int16_t adc = liftingsensor.getLastConversionResults();
    float voltage = liftingsensor.computeVolts(adc);
    return (84.22851f * voltage);
}

float readTailboardSensorRaw()
{
    int16_t adc = tailboardsensor.getLastConversionResults();
    float voltage = tailboardsensor.computeVolts(adc);
    return (61.0351f * voltage);
}

double calculate_alpha(double depthtarget, double tailboardangle)
{
    if (tailboardangle < 0.5)
        return 0.0;
    depthtarget = constrain(depthtarget, 40.0, 150.0);

    const double c0 = 14.8231841565, c1 = -2.5728483250, c2 = 2.5589109354, c3 = 0.0431399228,
                 c4 = 0.1178597710, c5 = -0.0222313485, c6 = 0.0032458064, c7 = -0.0068033308,
                 c8 = 0.0022842150, c9 = 0.0001798013;

    double a3 = c6;
    double a2 = c3 + c7 * tailboardangle;
    double a1 = c1 + c4 * tailboardangle + c8 * tailboardangle * tailboardangle;
    double a0 = c0 + c2 * tailboardangle + c5 * tailboardangle * tailboardangle + c9 * pow(tailboardangle, 3);

    double alpha = 20.0;
    for (int i = 0; i < 50; i++)
    {
        double f = a3 * pow(alpha, 3) + a2 * pow(alpha, 2) + a1 * alpha + (a0 - depthtarget);
        double df = 3.0 * a3 * pow(alpha, 2) + 2.0 * a2 * alpha + a1;
        if (fabs(df) < 1e-12)
            break;
        double delta = f / df;
        alpha = constrain(alpha - delta, 0.0, 40.0);
        if (fabs(delta) < 1e-8)
            break;
    }
    return alpha;
}

double calculate_estimate_depth(double liftingange, double tailboardangle)
{
    float estimate_depth = 0.0;
    if (tailboardangle < 0.5)
        return 0.0;
    estimate_depth = liftingange * 1.1 + tailboardangle * 2.0;
    return estimate_depth;
}

void driveActuator(float pwm_val)
{
    pwm_val = constrain(pwm_val, -PWM_MAX_ABS, PWM_MAX_ABS);
    if (pwm_val > 0.0f)
    {
        pwmIN1.writeScaled(pwm_val);
        pwmIN2.writeScaled(0.0f);
    }
    else if (pwm_val < 0.0f)
    {
        pwmIN1.writeScaled(0.0f);
        pwmIN2.writeScaled(abs(pwm_val));
    }
    else
    {
        pwmIN1.writeScaled(0.0f);
        pwmIN2.writeScaled(0.0f);
    }
}

void runController()
{
    float liftingangle = g_lifting_filtered;
    float tailangle = g_tail_filtered;
    float depth_target, Kp, Ki, Kd, K1, K2, K, tau1, L, alpha_manual;
    bool is_auto, dirty = false;

    if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(5)) == pdTRUE)
    {
        depth_target = g_depth_target;
        Kp = g_Kp;
        Ki = g_Ki;
        Kd = g_Kd;
        K1 = g_K1;
        K2 = g_K2;
        K = g_K;
        tau1 = g_tau1;
        L = g_L;
        is_auto = g_is_auto;
        alpha_manual = g_alpha_manual;
        dirty = g_model_dirty;
        g_model_dirty = false;
        xSemaphoreGive(g_mutex);
    }
    else
        return;

    float sp_raw = is_auto ? (float)calculate_alpha(depth_target, tailangle) : alpha_manual;
    g_sp_raw_val = sp_raw;

    ReferenceFilter::FilterOutput refOut = refFilter.update(sp_raw, DT);
    float alpha_ref = refOut.val, dot_alpha_ref = refOut.dot, ddot_alpha_ref = refOut.ddot;

    if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(5)) == pdTRUE)
    {
        g_setpoint = alpha_ref;
        g_dot_alpha_ref = dot_alpha_ref;   // Bổ sung dòng này
        g_ddot_alpha_ref = ddot_alpha_ref; // Bổ sung dòng này
        g_tailboardangle = tailangle;
        xSemaphoreGive(g_mutex);
    }

    if (dirty)
    {
        g_model_wd.K = -K; // u dương -> alpha giảm nên gain âm
        g_model_wd.tau1 = tau1;
        g_model_wd.L = L;
        g_model_wd.init(DT);
        g_model_nd.K = -K;
        g_model_nd.tau1 = tau1;
        g_model_nd.L = 0.0f;
        g_model_nd.init(DT);
    }

    static float u_prev = 0.0f;
    float y_m = g_model_nd.step(u_prev, DT);
    float y_m_d = g_model_wd.step(u_prev, DT);
    float y_pred = liftingangle + (y_m - y_m_d);

    static float y_pred_prev = 0.0f, e_int = 0.0f;
    float alpha_actual_dot = (y_pred - y_pred_prev) / DT;
    y_pred_prev = y_pred;
    float e = alpha_ref - y_pred, de = dot_alpha_ref - alpha_actual_dot;
    // e_int = constrain(e_int + e * DT, -500.0f, 500.0f);
    if (fabsf(e) > 0.5f) {
        e_int += e * DT;
    }
    e_int= constrain(e_int, -2.0f, 2.0f);

    float s = Kp * e + Ki * e_int + Kd * de;

    float safe_tau1 = (fabsf(tau1) < 1e-4f) ? 1e-4f : tau1;
    float a1 = 1.0f / safe_tau1, a2 = 0.0f, b = -K / safe_tau1;

    float u_eq = 0.0f;
    float safe_Kd = (fabsf(Kd) < 1e-4f) ? ((Kd >= 0.0f) ? 1e-4f : -1e-4f) : Kd;
    if (fabsf(b) > 1e-6f)
    {
        u_eq = (1.0f / b) * (a1 * alpha_actual_dot + a2 * y_pred + ddot_alpha_ref + (Kp / safe_Kd) * de + (Ki / safe_Kd) * e);
    }
    u_eq = constrain(u_eq, -PWM_MAX_ABS, PWM_MAX_ABS);

    static float u_sw_int = 0.0f;
    float sign_s = (s > 0.0f) ? 1.0f : ((s < 0.0f) ? -1.0f : 0.0f);

    // Super-Twisting Algorithm cho hệ gain âm (b < 0):
    // u_sw_int tích lũy cùng dấu với s để khi chia cho b_Kd (âm) sẽ ra u cùng chiều cần thiết
    u_sw_int = constrain(u_sw_int + K2 * sign_s * DT, -2.0f, 2.0f);

    float u_sw = 0.0f;
    float b_Kd = b * safe_Kd;
    if (fabsf(b_Kd) > 1e-6f)
    {
        // Loại bỏ dấu trừ trước K1: u_sw = (1/b_Kd) * (K1*sqrt|s|*sign_s + u_sw_int)
        u_sw = (1.0f / b_Kd) * (K1 * sqrtf(fabsf(s)) * sign_s + u_sw_int);
    }
    u_sw = constrain(u_sw, -PWM_MAX_ABS, PWM_MAX_ABS);

    float u_out = constrain(u_eq + u_sw, -PWM_MAX_ABS, PWM_MAX_ABS);
    driveActuator(u_out);
    u_prev = u_out;
    float estimatedepth = (float)calculate_estimate_depth(liftingangle, tailangle);
    if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(5)) == pdTRUE)
    {
        g_liftingangle = liftingangle;
        g_e = e;
        g_de = de;
        g_s = s;
        g_u = u_out;
        g_u_eq = u_eq;
        g_u_sw = u_sw;
        g_u_sw_i = u_sw_int;
        g_estimate_depth = estimatedepth;
        g_dot_alpha_ref = dot_alpha_ref;
        g_ddot_alpha_ref = ddot_alpha_ref;
        xSemaphoreGive(g_mutex);
    }
}

void sensorReadTask(void *param)
{
    TickType_t xLastWake = xTaskGetTickCount();
    while (true)
    {
        float lifting_raw = g_lift_offset - readLiftingSensorRaw();
        float tail_raw = g_tail_offset - readTailboardSensorRaw();
        g_lifting_filtered = filterLifting.update(lifting_raw);
        g_tail_filtered = filterTailboard.update(tail_raw);
        g_lifting_raw_val = lifting_raw;
        g_tail_raw_val = tail_raw;
        vTaskDelayUntil(&xLastWake, SENSOR_DT_TICKS);
    }
}

void controlTask(void *param)
{
    TickType_t xLastWake = xTaskGetTickCount();
    uint32_t count = 0;
    uint32_t last_report_ms = millis();

    while (true)
    {
        runController();
        count++;

        uint32_t now = millis();
        if (now - last_report_ms >= 1000)
        {
            float freq = (float)count * 1000.0f / (float)(now - last_report_ms);
            Serial.printf("[DEBUG] controlTask Frequency: %.2f Hz\n", freq);
            count = 0;
            last_report_ms = now;
        }

        vTaskDelayUntil(&xLastWake, DT_TICKS);
    }
}

#ifdef MONITOR_WIFI

struct __attribute__((packed)) SettingsPayload
{
    float depth_target; // Độ sâu mục tiêu (mm)
    float Kp;
    float Ki;
    float Kd;
    float K1;
    float K2;
    float K;
    float tau1;
    float tau2;
    float L;
    float lift_offset;
    float tail_offset;
    float is_auto;      // 0=Manual, 1=Auto
    float alpha_manual; // Alpha target thủ công
};

struct __attribute__((packed)) ControllerStatus
{
    float liftingangle;
    float tailboardangle;
    float setpoint;
    float depth_target;
    float e;
    float de;
    float s;
    float u;
    float estimate_depth;
    float is_auto;
};

void broadcastData()
{
    if (ws.count() == 0)
        return;

    ControllerStatus status;
    if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(5)) != pdTRUE)
        return;
    status.liftingangle = g_liftingangle;
    status.tailboardangle = g_tailboardangle;
    status.setpoint = g_setpoint;
    status.depth_target = g_depth_target;
    status.e = g_e;
    status.de = g_de;
    status.s = g_s;
    status.u = g_u;
    status.is_auto = g_is_auto ? 1.0f : 0.0f;
    xSemaphoreGive(g_mutex);

    ws.binaryAll((uint8_t *)&status, sizeof(ControllerStatus));
}

void webBroadcastTask(void *param)
{
    TickType_t xLastWake = xTaskGetTickCount();
    while (true)
    {
        broadcastData();
        ws.cleanupClients();
        vTaskDelayUntil(&xLastWake, WEB_DT_TICKS);
    }
}

void onWsEvent(AsyncWebSocket *server, AsyncWebSocketClient *client, AwsEventType type, void *arg, uint8_t *data, size_t len)
{
    if (type == WS_EVT_DATA)
    {
        AwsFrameInfo *info = (AwsFrameInfo *)arg;
        if (info->final && info->index == 0 && info->len == len && info->opcode == WS_BINARY)
        {
            if (len == sizeof(SettingsPayload))
            {
                SettingsPayload *p = (SettingsPayload *)data;
                if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(10)) == pdTRUE)
                {
                    g_depth_target = constrain(p->depth_target, 40.0f, 150.0f);
                    g_Kp = constrain(p->Kp, 0.01f, 200.0f);
                    g_Ki = constrain(p->Ki, 0.0f, 100.0f);
                    g_Kd = constrain(p->Kd, 0.001f, 100.0f);
                    g_K1 = constrain(p->K1, 0.01f, 500.0f);
                    g_K2 = constrain(p->K2, 0.01f, 500.0f);
                    g_lift_offset = constrain(p->lift_offset, 0.0f, 345.0f);
                    g_tail_offset = constrain(p->tail_offset, 0.0f, 250.f);
                    g_is_auto = (p->is_auto > 0.5f);
                    g_alpha_manual = constrain(p->alpha_manual, 0.0f, 40.0f);

                    if (p->K != g_K || p->tau1 != g_tau1 || p->L != g_L)
                    {
                        g_K = constrain(p->K, 0.01f, 100.0f);
                        g_tau1 = constrain(p->tau1, 0.01f, 100.0f);
                        g_L = constrain(p->L, 0.0f, 10.0f);
                        g_model_dirty = true;
                    }
                    xSemaphoreGive(g_mutex);
                    Serial.println("[WS] Binary Settings Updated");
                }
            }
        }
        else if (info->opcode == WS_TEXT)
        {
            // Vẫn giữ lại xử lý text JSON nếu cần test thủ công
        }
    }
}

void setupRoutes()
{
    ws.onEvent(onWsEvent);
    server.addHandler(&ws);

    // Captive Portal: Điều hướng các yêu cầu phổ biến của Android/iOS
    server.on("/generate_204", HTTP_GET, [](AsyncWebServerRequest *request)
              { request->redirect("/"); }); // Android
    server.on("/fwlink", HTTP_GET, [](AsyncWebServerRequest *request)
              { request->redirect("/"); }); // Windows

    server.on("/", HTTP_GET, [](AsyncWebServerRequest *req)
              {
        AsyncWebServerResponse *response = req->beginResponse(LittleFS, "/index.html", "text/html");
        req->send(response); });

    server.onNotFound([](AsyncWebServerRequest *req)
                      {
        // Nếu yêu cầu không khớp với bất kỳ route nào, điều hướng về trang chủ (cần thiết cho Captive Portal)
        req->redirect("/"); });
}
#endif

// ────────────────────────────────────────────────────────────────────
//  SERIAL TELEMETRY & TUNING
// ────────────────────────────────────────────────────────────────────
volatile float g_fc_lifting = 25.0f, g_fc_tailboard = 25.0f, g_omega_ref = 10.0f;

struct __attribute__((packed)) SerialTelemetry
{
    uint16_t header = 0xAA55;
    float values[32]; // 0-5: alpha_r/f, th_r/f, sp_r/f | 6-10: dtar, e, de, s, u | 11-15: Kp, Ki, Kd, K1, K2 | 16-18: K, t1, L | 19-21: fc_L, fc_T, om_R | 22-23: L_off, T_off | 24: is_auto | 25: alpha_man | 26: u_eq | 27: u_sw | 28: u_sw_i | 29: est_depth | 30: dot_ref | 31: ddot_ref
    uint8_t checksum;
};

struct __attribute__((packed)) SerialCommand
{
    uint16_t header;
    float values[16]; // dtar, Kp, Ki, Kd, K1, K2, K, t1, L, fc_L, fc_T, om_R, L_off, T_off, is_auto, alpha_man
    uint8_t checksum;
};

#ifdef DEBUG_SERIAL
void serialTuningTask(void *param)
{
    TickType_t xLastWake = xTaskGetTickCount();
    while (true)
    {
        SerialTelemetry msg;
        if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(5)) == pdTRUE)
        {
            float *v = msg.values;
            v[0] = g_lifting_raw_val;
            v[1] = g_liftingangle;
            v[2] = g_tail_raw_val;
            v[3] = g_tailboardangle;
            v[4] = g_sp_raw_val;
            v[5] = g_setpoint;
            v[6] = g_depth_target;
            v[7] = g_e;
            v[8] = g_de;
            v[9] = g_s;
            v[10] = g_u;
            v[11] = g_Kp;
            v[12] = g_Ki;
            v[13] = g_Kd;
            v[14] = g_K1;
            v[15] = g_K2;
            v[16] = g_K;
            v[17] = g_tau1;
            v[18] = g_L;
            v[19] = g_fc_lifting;
            v[20] = g_fc_tailboard;
            v[21] = g_omega_ref;
            v[22] = g_lift_offset;
            v[23] = g_tail_offset;
            v[24] = g_is_auto ? 1.0f : 0.0f;
            v[25] = g_alpha_manual;
            v[26] = g_u_eq;
            v[27] = g_u_sw;
            v[28] = g_u_sw_i;
            v[29] = g_estimate_depth;
            v[30] = g_dot_alpha_ref;
            v[31] = g_ddot_alpha_ref;
            xSemaphoreGive(g_mutex);
        }
        msg.checksum = calculateChecksum((uint8_t *)&msg.values, sizeof(msg.values));
        Serial.write((uint8_t *)&msg, sizeof(msg));

        if (Serial.available() >= sizeof(SerialCommand))
        {
            uint8_t buf[sizeof(SerialCommand)];
            if (Serial.peek() == 0x55)
            {
                Serial.readBytes(buf, sizeof(SerialCommand));
                SerialCommand *cmd = (SerialCommand *)buf;
                if (cmd->header == 0xAA55 && calculateChecksum((uint8_t *)&cmd->values, sizeof(cmd->values)) == cmd->checksum)
                {
                    if (xSemaphoreTake(g_mutex, pdMS_TO_TICKS(10)) == pdTRUE)
                    {
                        float *cv = cmd->values;
                        g_depth_target = cv[0];
                        g_Kp = cv[1];
                        g_Ki = cv[2];
                        g_Kd = cv[3];
                        g_K1 = cv[4];
                        g_K2 = cv[5];
                        g_K = cv[6];
                        g_tau1 = cv[7];
                        g_L = cv[8];
                        if (g_fc_lifting != cv[9] || g_fc_tailboard != cv[10] || g_omega_ref != cv[11])
                        {
                            g_fc_lifting = cv[9];
                            g_fc_tailboard = cv[10];
                            g_omega_ref = cv[11];
                            filterLifting.init(g_fc_lifting, 500.0f);
                            filterTailboard.init(g_fc_tailboard, 500.0f);
                            refFilter.init(g_omega_ref, 1.0f);
                        }
                        g_lift_offset = cv[12];
                        g_tail_offset = cv[13];
                        g_is_auto = (cv[14] > 0.5f);
                        g_alpha_manual = cv[15];
                        g_model_dirty = true;
                        xSemaphoreGive(g_mutex);
                    }
                }
            }
            else
                Serial.read();
        }
        vTaskDelayUntil(&xLastWake, SERIALDEBUG_DT_TICKS);
    }
}
#endif

void setup()
{
    setCpuFrequencyMhz(240);
    delay(500);
    Serial.begin(115200);
    delay(500);
    Serial.println("\n=== Tractor SOIPDT + Super-Twisting SMC ===");
    Wire.begin(PIN_SDA, PIN_SCL);
    Wire.setClock(400000);
    if (liftingsensor.begin(0x48))
    {
        liftingsensor.setGain(GAIN_ONE);
        liftingsensor.setDataRate(RATE_ADS1115_860SPS);
        liftingsensor.startADCReading(ADS1X15_REG_CONFIG_MUX_DIFF_0_1, true);
    }
    if (tailboardsensor.begin(0x49))
    {
        tailboardsensor.setGain(GAIN_ONE);
        tailboardsensor.setDataRate(RATE_ADS1115_860SPS);
        tailboardsensor.startADCReading(ADS1X15_REG_CONFIG_MUX_DIFF_0_1, true);
    }
    ESP32PWM::allocateTimer(0);
    ESP32PWM::allocateTimer(1);
    ESP32PWM::allocateTimer(2);
    ESP32PWM::allocateTimer(3);
    pwmIN1.attachPin(PIN_MOTOR_IN1, PWM_FREQ_HZ, PWM_BITS);
    pwmIN2.attachPin(PIN_MOTOR_IN2, PWM_FREQ_HZ, PWM_BITS);
    g_mutex = xSemaphoreCreateMutex();
    filterLifting.init(25.0f, 500.0f);
    filterTailboard.init(25.0f, 500.0f);
    refFilter.init(10.0f, 1.0f);
    g_model_wd.K = g_K;
    g_model_wd.tau1 = g_tau1;
    g_model_wd.L = g_L;
    g_model_wd.init(DT);
    g_model_nd.K = g_K;
    g_model_nd.tau1 = g_tau1;
    g_model_nd.L = 0.0f;
    g_model_nd.init(DT);
    g_model_dirty = false;

#ifdef MONITOR_WIFI

    if (LittleFS.begin(true))
    {
        Serial.println("[FS] LittleFS mounted successfully");
    }
    else
    {
        Serial.println("[FS] LittleFS mount failed");
    }

    // WiFi AP
    WiFi.mode(WIFI_AP);
    WiFi.softAP(AP_SSID, AP_PASS);
    dnsServer.start(53, "*", WiFi.softAPIP());
    Serial.printf("[WiFi]  SSID: %-20s  IP: %s\n",
                  AP_SSID, WiFi.softAPIP().toString().c_str());

    // Web server
    setupRoutes();
    server.begin();
    Serial.println("[HTTP]  Server OK – port 80");
    // xTaskCreatePinnedToCore(wsStatusTask, "WSStatus", 4096, nullptr, 1, nullptr, 0);
    xTaskCreatePinnedToCore(webBroadcastTask, "Web_WS", 4096, nullptr, 1, nullptr, 0);
#endif

    xTaskCreatePinnedToCore(sensorReadTask, "SensRead", 4096, nullptr, 3, nullptr, 1);
    xTaskCreatePinnedToCore(controlTask, "SMC_Ctrl", 8192, nullptr, 2, nullptr, 1);
#ifdef DEBUG_SERIAL
    xTaskCreatePinnedToCore(serialTuningTask, "SerTune", 4096, nullptr, 1, nullptr, 0);
#endif
    Serial.println("=== Ready! Serial Tuning Active ===\n");
}

void loop() { vTaskDelay(pdMS_TO_TICKS(100)); }
