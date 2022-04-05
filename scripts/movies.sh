#!/bin/bash

module load ffmpeg

frame_rate='15'
video_name='OP_MWLMC5_b0_over_density.mp4'
file_name='../mwlmc5_over_density/MWLMC5_over_density_100M_b0_OM3_'
ffmpeg  -r $frame_rate -i $file_name%3d.png -s 1920x1080  -crf 15 -pix_fmt yuv420p $video_name

