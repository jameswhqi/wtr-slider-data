#!/usr/bin/env python3

import sys
import json
import os
import os.path

data = json.loads(
    sys.stdin.buffer.read(int(os.environ["CONTENT_LENGTH"])).decode("utf-8")
)
inCgiBin = os.path.basename(os.getcwd()) == "cgi-bin"
with open(
    f"{'../' if inCgiBin else ''}data/{data['client']['workerId']}.json", "w"
) as f:
    json.dump(data, f, indent=2)

print("Content-Type:text/plain\n\n")
print("Success")
