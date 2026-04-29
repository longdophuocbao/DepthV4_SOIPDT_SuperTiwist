import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    df = pd.read_csv('tractor_log_20260428_135513.csv')
    
    # Calculate errors
    e = df['e']
    sp = df['sp_ref']
    y = df['y_pred']
    u = df['u']
    
    rmse = np.sqrt(np.mean(e**2))
    max_err = np.max(np.abs(e))
    mae = np.mean(np.abs(e))
    print("=== ERROR METRICS ===")
    print(f"RMSE: {rmse:.4f}")
    print(f"Max Error: {max_err:.4f}")
    print(f"MAE: {mae:.4f}")
    
    # Control effort
    u_rms = np.sqrt(np.mean(u**2))
    u_max = np.max(np.abs(u))
    print("\n=== CONTROL EFFORT ===")
    print(f"Control RMS: {u_rms:.4f}")
    print(f"Control Max: {u_max:.4f}")
    
    # Sliding surface
    s = df['s']
    s_rms = np.sqrt(np.mean(s**2))
    s_max = np.max(np.abs(s))
    print("\n=== SLIDING SURFACE ===")
    print(f"S RMS: {s_rms:.4f}")
    print(f"S Max: {s_max:.4f}")
    
    print("\n=== CURRENT PARAMS ===")
    for col in ['Kp', 'Ki', 'Kd', 'K1', 'K2', 'K', 'tau1', 'L']:
        if col in df.columns:
            print(f"{col}: {df[col].iloc[-1]}")
            
    # Try to plot and save
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1)
    plt.plot(y.values, label='y_pred')
    plt.plot(sp.values, label='sp_ref', linestyle='--')
    plt.legend()
    plt.title('Tracking Performance')
    
    plt.subplot(3, 1, 2)
    plt.plot(e.values, label='Error (e)')
    plt.plot(s.values, label='Sliding Surface (s)', alpha=0.7)
    plt.legend()
    plt.title('Error and Sliding Surface')
    
    plt.subplot(3, 1, 3)
    plt.plot(u.values, label='Control Signal (u)')
    plt.plot(df['u_eq'].values, label='u_eq', alpha=0.7)
    plt.plot(df['u_sw'].values, label='u_sw', alpha=0.7)
    plt.legend()
    plt.title('Control Effort')
    
    plt.tight_layout()
    plt.savefig('plot.png')
    print("\nPlot saved as plot.png")
except Exception as ex:
    print(f"Error: {ex}")
