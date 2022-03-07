from math import cos, sin

def rotate_by_angle(x, y, angle):
    x2 = x*cos(angle) - y*sin(angle)
    y2 = x*sin(angle) + y*cos(angle)
    return [x2, y2]

