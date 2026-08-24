# %%
import os
import json

folder = os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026', 'processed')

for f, _, fns in os.walk(folder):
    for fn in fns:
        if fn=='protocol.json':
            print('processing %s' % f)
            with open(os.path.join(f, fn), 'r') as r:
                protocol = r.read()

            if len(protocol.split("\"Protocol-2\": \"contrast-sensitivity.json\""))>0:
                print('BEFORE:')
                print(protocol[:600])
                print()

                protocol = protocol.replace(
                    "\"Protocol-2\": \"contrast-sensitivity.json\"",
                    "\"Protocol-2\": \"tuning-low-contrast.json\""
                )
                protocol = protocol.replace(
                    "\"Protocol-3\": \"tuning-low-contrast.json\"",
                    "\"Protocol-3\": \"tuning-mid-contrast.json\""
                )
                protocol = protocol.replace(
                    "\"Protocol-4\": \"tuning-mid-contrast.json\"",
                    "\"Protocol-4\": \"tuning-high-contrast.json\""
                )
                protocol = protocol.replace(
                    "\"Protocol-5\": \"tuning-high-contrast.json\"",
                    "\"Protocol-5\": \"contrast-sensitivity.json\""
                )
                print('AFTER:')
                print(protocol[:600])
                print()
                # over-write
                with open(os.path.join(f, fn), 'w') as r:
                    r.write(protocol)



# %%
