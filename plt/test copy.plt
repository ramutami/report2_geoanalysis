system("rm -r ../imageout/*")

# fortranファイルの出力を見て、ここを書き換えること！
    first_filenum = 1
    last_filenum = 730
#

set terminal pngcairo size 800,600
set output "../imageout/vect.png"

set xlabel "経度[deg]"
set xrange [0:60]
set yrange [20:40]
set ylabel "緯度[deg]"

set cblabel "|U| [m/s]"
set colorbox
set palette defined ( \
    0.00 "blue", \
    0.05 "cyan", \
    0.10 "green", \
    0.20 "yellow", \
    0.30 "orange", \
    0.50 "red" \
)
set cbrange [0:0.25]

set size ratio -2
scale = 10.0

mag(u,v) = sqrt((u)**2 + (v)**2)
eps = 1e-4

do for [n=first_filenum:last_filenum]{
    imageout_name = sprintf("../imageout/vect%06d.png", n)
    file_in_name = sprintf("../dataout/dataout%06d.dat", n)
    set output imageout_name
    time_dat = real(system(sprintf("awk 'NR==1{print $3}' %s", file_in_name)))
    step_num = int(system(sprintf("awk 'NR==2{print $4}' %s", file_in_name)))
    set title sprintf("t = %.1f day, step = %d", time_dat, step_num)

    plot file_in_name using 1:2:(mag($4,$5) > eps ? $4*scale : 1/0):(mag($4,$5) > eps ? $5*scale : 1/0):(mag($4,$5)) with vectors head filled size screen 0.004,8,25 lc palette lw 0.4 notitle
    print sprintf("plotting %d / %d", n, last_filenum)
}


system("rm ../videoout/vect.mp4")
system("ffmpeg -framerate 30 \
-i ../imageout/vect%06d.png \
-c:v libx264 \
-pix_fmt yuv420p \
-movflags +faststart \
../videoout/vect.mp4")


#system("rm -r ../imageout/*")
#system("rm -r ../dataout/*")