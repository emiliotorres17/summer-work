#!/usr/bin/env python3
import os
import sys
from subprocess import call


if __name__ == '__main__':
    path = 'media'
    for filename in os.listdir(path):
        if (filename.endswith(".png")): #or .avi, .mpeg, whatever.
            call(['ffmpeg', '-i', 'media/plot-%d0.png', 'output.mp4'])
        else:
            continue
