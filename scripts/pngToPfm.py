# png_to_pfm.py
import numpy as np
import imageio.v3 as iio

def write_pfm(filename, image):
    height, width, _ = image.shape
    with open(filename, 'wb') as f:
        f.write(b'PF\n')
        f.write(f"{width} {height}\n".encode())
        f.write(b'-1.0\n')  # Little-endian
        f.write(image[::-1].astype(np.float32).tobytes())

img = iio.imread('input.png').astype(np.float32) / 255.0  # Normalize to [0, 1]
write_pfm('input.pfm', img)
