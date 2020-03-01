#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

values_train = []
for line in open(sys.argv[1]):
    fields=line.strip('\r\n').split()
    for value in fields:
        values_train.append(float(value))

values_test = []
for line in open(sys.argv[2]):
    fields=line.strip('\r\n').split()
    for value in fields:
        values_test.append(float(value))

values_ref = []
for line in open(sys.argv[3]):
    fields=line.strip('\r\n').split()
    for value in fields:
        values_ref.append(float(value))

fig, (ax1, ax2, ax3) = plt.subplots(1,3)
fig.suptitle('log TPM')
ax1.set_title('Train')
ax1.hist(values_train, bins=100)
ax1.set_xlim(-10, 10)
ax2.set_title('Test')
ax2.hist(values_test, bins=100)
ax2.set_xlim(-10, 10)
ax3.set_title('Ref')
ax3.hist(values_ref, bins=100)
ax3.set_xlim(-10, 10)
fig.savefig(sys.argv[4])
plt.close(fig)
