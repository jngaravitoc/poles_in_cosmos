#!/bin/bash

frame_rate='15'
video_name='OP_MWLMC5_b0_dens_0.mp4'
ffmpeg  -r $frame_rate -i ./OP_MWLMC5_100M_b0_OM3_%3d.png -s 1920x1080 -vcodec libx264 -crf 15 -pix_fmt yuv420p $video_name

