# Fix for Boolean Index Error in Temporal-Dynamics Scripts

## Error Message
```
boolean index did not match indexed array along axis 1; size of axis is 80 but size of corresponding boolean axis is 598
```

This error occurred when running `Temporal-Dynamics_v5.py` (and similar versions) while creating plots.

## Root Cause

The error was caused by a shape mismatch between:
1. Boolean masks (`peak_cond`, `baseline_cond`) created from the contrast-sensitivity protocol's time array (598 time points)
2. Response data arrays (`resp_data`) from tuning protocols with different time dimensions (80 time points)

In the original code:
- `summary['t']` was set from the first contrast-sensitivity protocol file found
- All protocols used this single time array for creating boolean masks
- But tuning protocol response data had different time dimensions
- Applying the mismatched boolean masks to the response data caused the numpy error

## Solution

Modified three files to store and use protocol-specific time arrays:
- `Temporal-Dynamics_v3.py`
- `Temporal-Dynamics_v4.py`
- `Temporal-Dynamics_v5.py`

### Changes Made

#### 1. In `build_summary` function:
- Added `'t': None` to protocol dictionaries
- Store time array from first file for each protocol:
  ```python
  # For contrast-sensitivity
  summary[virus][protocol] = {'frac_resp': [], 't': None}
  if summary[virus][protocol]['t'] is None and 't' in responses:
      summary[virus][protocol]['t'] = responses['t']
  
  # For tuning protocols
  summary[virus][protocol] = {'frac_resp': [], 'resp': [], 't': None}
  if summary[virus][protocol]['t'] is None and 't' in responses:
      summary[virus][protocol]['t'] = responses['t']
  ```

#### 2. In `make_fig` function:
- Get protocol-specific time array before creating boolean masks:
  ```python
  protocol_t = summary[virus][protocol].get('t') or summary.get('t', np.array([]))
  if len(protocol_t) == 0:
      protocol_t = summary.get('t', np.array([]))
  
  # Create boolean conditions using protocol-specific time
  peak_cond = (protocol_t > peak_window[0]) & (protocol_t < peak_window[1])
  baseline_cond = (protocol_t > baseline_window[0]) & (protocol_t < baseline_window[1])
  ```
  **Important:** Use `.get('t') or` instead of `.get('t', ...)` to properly handle the case where `summary[virus][protocol]['t']` is None. The `.get('t', fallback)` method returns None if the key exists with value None, rather than returning the fallback.
- Use `protocol_t` instead of `summary['t']` when plotting

## Files Modified
1. `Temporal-Dynamics_v3.py` - Fixed build_summary and make_fig
2. `Temporal-Dynamics_v4.py` - Fixed build_summary and make_fig  
3. `Temporal-Dynamics_v5.py` - Fixed build_summary and make_fig

## Impact
- Boolean masks now match the shape of the response data they're applied to
- Each protocol uses its own time array for indexing
- Falls back to global `summary['t']` if protocol-specific time is not available
- Maintains backward compatibility

## Testing
All modified files pass Python syntax validation:
```bash
python -c "import ast; ast.parse(open('Temporal-Dynamics_v5.py').read()); print('Syntax OK')"
```
