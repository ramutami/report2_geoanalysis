system("rm -r ../imageout/*")

# fortranファイルの出力を見て、ここを書き換えること！
    first_filenum = 1
    last_filenum = 730
#

set terminal pngcairo size 800,600
set output "../imageout/heatmap.png"

set xlabel "経度[deg]"
set xrange [0:60]
set yrange [20:40]
set ylabel "緯度[deg]"
set cblabel "h[m]"
set cbrange [-2:10]

set view map
set size ratio -2
set pm3d map
set colorbox

set palette defined ( \
    -2 "purple",\
    0 "white", \
    2 "blue", \
    4 "green", \
    6 "yellow", \
    8 "orange", \
    10 "red" \
)

do for [n=first_filenum:last_filenum]{
    imageout_name = sprintf("../imageout/h%06d.png", n)
    file_in_name = sprintf("../dataout/dataout%06d.dat", n)
    set output imageout_name
    time_dat = real(system(sprintf("awk 'NR==1{print $3}' %s", file_in_name)))
    step_num = int(system(sprintf("awk 'NR==2{print $4}' %s", file_in_name)))
    set title sprintf("t = %.1f day, step = %d", time_dat, step_num)

    splot file_in_name using 1:2:3 with pm3d notitle
    print sprintf("plotting %d / %d", n, last_filenum)
}


system("rm ../videoout/h.mp4")
system("ffmpeg -framerate 30 \
-i ../imageout/h%06d.png \
-c:v libx264 \
-pix_fmt yuv420p \
-movflags +faststart \
../videoout/h.mp4")
#system("rm -r ../imageout/*")

#system("rm -r ../dataout/*")