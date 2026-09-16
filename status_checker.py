import matplotlib.pyplot as plt
import numpy as np

# Load the file
data = np.loadtxt("status_checker,N=100.dat")

x = data[:, 0]
y = data[:, 1]
z = data[:, 2]
status = data[:, 3]

# Create boolean masks for status filtering
unoccupied = status == 0
occupied = (status == 1) | (status == 2)
travelling = status == 3

fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(projection="3d")
ax = fig.add_subplot()

# Plot unoccupied points (Green, lower opacity)
ax.scatter(
    x[unoccupied],
    y[unoccupied],
    # z[unoccupied],
    c="green",
    alpha=0.05,
    marker="o",
    label="Unoccupied",
)

# Plot occupied & travelling points (Red, higher opacity)
ax.scatter(
    x[occupied],
    y[occupied],
    # z[occupied],
    c="red",
    alpha=0.8,
    marker="o",
    label="Occupied / Travelling",
)


ax.scatter(
    x[occupied],
    y[occupied],
    # z[occupied],
    c="blue",
    alpha=0.10,
    marker="o",
    label="Targets",
)

ax.set_xlabel("X")
ax.set_ylabel("Y")
# ax.set_zlabel("Z")
ax.legend()
plt.show()