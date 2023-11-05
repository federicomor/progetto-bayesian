# sudo apt install ffmpeg
if [ -z "$1" ]; then
	echo "Error: provide images-folder as argument (ie call 'bash create_video.sh FOLDER')."
	exit 1
fi

cd "$1" || exit 1
ffmpeg -framerate 5 -pattern_type glob -i '*.jpeg' -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4
