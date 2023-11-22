# sudo apt install imagemagick
if [ -z "$1" ]; then
	echo "Error: provide images-folder as argument (ie call 'bash create_gif.sh FOLDER')."
	exit 1
fi

cd "$1" || exit 1
convert -resize 90% -delay 15 -verbose -loop 0 *.jpeg out.gif
