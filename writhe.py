import numpy as np


def read_3col(filename, num_bp, num_steps):
    """
    Reads a 3col coordinate file and splits it by timestep
    """
    coords = np.loadtxt(filename)
    x = []
    t = 0
    for i in range(num_steps):
        x.append(coords[t*num_bp:(t+1)*num_bp, :])
        t += 1
    return np.array(x)


def writhe(x, t, l, axis=2):
    """
    Calculates write for a single timestep
    """
    shape = np.shape(x[t])
    y = np.zeros((shape[0]+1, shape[1]))
    # Read one bp-step
    y[:shape[0], :] = x[t]
    # Add the head at the bottom
    y[shape[0], :] = x[t, 0]
    result = 0
    for j in range(l):
        for k in range(j):
            tangents_j = y[j+1] - y[j]
            tangents_k = y[k+1] - y[k]
            # Vector joining j and k
            vector = y[j] - y[k]
            # Discretised Gauss integral
            # Add each individual contribution from a pair of tangent vectors
            result += (np.dot(vector, np.cross(tangents_j, tangents_k)) /
                       (np.linalg.norm(vector)**3 * 2 * np.pi))
    return result


def main(name, num_bp, num_steps):
    # Read file
    coords = read_3col(name + '/C1.3col', num_bp, num_steps)
    length = len(coords[0])
    # Calculate writhe for num_steps timesteps
    wr = []
    for t in range(num_steps):
        print(f"\r\tStep {t}...", end=" ")
        wr.append([t+1, writhe(coords, t, length)])
    wr = np.array(wr)
    np.savetxt(name+'/writhe.ser', wr, fmt='%5d %9.4f')
    print("Done!")
    return wr
