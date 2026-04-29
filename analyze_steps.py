import pandas as pd
import numpy as np

df = pd.read_csv('tractor_log_20260428_135513.csv')

df['sp_raw_diff'] = df['alpha_manual'].diff()
step_indices = df.index[np.abs(df['sp_raw_diff']) > 0.5].tolist()

print(f"Total rows: {len(df)}")
print(f"Step changes in alpha_manual at indices: {step_indices}")

for idx in step_indices:
    old_sp = df['alpha_manual'].iloc[idx-1]
    new_sp = df['alpha_manual'].iloc[idx]
    print(f"\n--- Step from {old_sp} to {new_sp} ---")
    
    # Evaluate response over the next 1500 samples (~3 seconds)
    end_idx = min(idx + 1500, len(df))
    segment = df.iloc[idx:end_idx]
    
    # e is y_pred - sp_ref
    # tracking error relative to target sp is y_pred - alpha_manual
    # tracking error relative to smoothed ref is e
    e_smooth = segment['e']
    e_target = segment['y_pred'] - segment['alpha_manual']
    
    max_err_smooth = np.max(np.abs(e_smooth))
    rmse_smooth = np.sqrt(np.mean(e_smooth**2))
    
    print(f"Tracking vs smoothed ref -> Max Error: {max_err_smooth:.4f}, RMSE: {rmse_smooth:.4f}")
    
    # Calculate overshoot based on actual target
    y_pred = segment['y_pred']
    if new_sp > old_sp:
        overshoot = np.max(y_pred) - new_sp
    else:
        overshoot = new_sp - np.min(y_pred)
    overshoot_percent = max(0, overshoot / np.abs(new_sp - old_sp) * 100)
    print(f"Overshoot: {overshoot:.4f} ({overshoot_percent:.2f}%)")
    
    # Settling time to target
    threshold = 0.5
    settled_indices = segment.index[np.abs(e_target) < threshold].tolist()
    if len(settled_indices) > 0:
        outside = segment.index[np.abs(e_target) >= threshold].tolist()
        if len(outside) > 0:
            settle_idx = outside[-1] + 1
            if settle_idx < end_idx:
                settle_time = (settle_idx - idx) * 0.002
                print(f"Settling time to target (+/- {threshold}): {settle_time:.4f} s")
            else:
                print("Did not settle within window")
        else:
            print("Always within threshold")
    else:
        print("Did not reach threshold")
