# sudo apt install ffmpeg
ffmpeg -framerate 5 -pattern_type glob -i '*.jpeg' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4