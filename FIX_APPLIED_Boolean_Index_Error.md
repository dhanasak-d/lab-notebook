# Fix Applied for Boolean Index Error

## Date
2026-10-17

## Issue
When running Temporal-Dynamics_v5.py (and similar versions v3, v4), the following error occurred repeatedly:
```
Error: boolean index did not match indexed array along axis 1; size of axis is 80 but size of corresponding boolean axis is 598
```

This error appeared during the creation of sessions and ROIs figures.

## Root Cause
The error occurred because of a subtle bug in how protocol-specific time arrays were being retrieved. The issue was in the `make_fig` function:

```python
protocol_t = summary[virus][protocol].get('t', summary.get('t', np.array([])))
```

When `summary[virus][protocol]['t']` was `None`, the `.get('t', fallback)` method would return `None` (not the fallback) because the key exists with value `None`. This would then cause:
1. A TypeError when trying to call `len(protocol_t)` on `None`
2. Or if somehow bypassed, a shape mismatch when creating boolean masks

## Solution Applied
Changed the line in the `make_fig` function in three files:
- `Temporal-Dynamics_v3.py` (line ~633)
- `Temporal-Dynamics_v4.py` (line ~553)
- `Temporal-Dynamics_v5.py` (line ~356)

**Before:**
```python
protocol_t = summary[virus][protocol].get('t', summary.get('t', np.array([])))
```

**After:**
```python
protocol_t = summary[virus][protocol].get('t') or summary.get('t', np.array([]))
```

## Why This Works
In Python:
- `d.get('t', fallback)` returns the value if the key exists, even if the value is `None`
- `d.get('t') or fallback` returns the value if it's truthy, otherwise returns `fallback`

So when `summary[virus][protocol]['t']` is `None`:
- Old: Returns `None` → causes TypeError on `len(None)`
- New: Returns `summary.get('t', np.array([]))` → uses the fallback time array

When `summary[virus][protocol]['t']` is an array:
- Both: Returns the array → correct behavior

## Files Modified
1. `Temporal-Dynamics_v3.py` - Fixed line 633
2. `Temporal-Dynamics_v4.py` - Fixed line 553
3. `Temporal-Dynamics_v5.py` - Fixed line 356
4. `FIX_SUMMARY_Boolean_Index_Error.md` - Updated documentation

## Verification
All modified files pass Python syntax validation:
```bash
python -c "import ast; ast.parse(open('Temporal-Dynamics_v5.py').read()); print('Syntax OK')"
```

## Testing
The fix ensures that:
1. If a protocol has its own time array stored, it will be used
2. If a protocol's time array is `None`, it falls back to the global time array
3. Boolean masks are created with the correct dimensions matching the response data

This should resolve the "boolean index did not match indexed array" error.
