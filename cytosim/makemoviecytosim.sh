# Neha, this generates a movie from the ppm files and then deletes all the ppm files.
# Laste modified, 15 April 2017




## For uncompressed version
#  ffmpeg -f image2 -i image%04d.ppm -an -vcodec rawvideo -y output.avi

# to rEduces typical 270 MB avi movie to 5 MB movie!! :)

cellRad=$1
mtdyn=$2
mottype=$3

#dyn=$3
#kin=$4
#$lindensity=$5



#filename=R$cellRad'_'MtDyn$mtdyn'_'noMotors
#filename=R$cellRad'_'MtDyn$mtdyn'_'cortical_dyn_diff_kin


#filename=R$cellRad'_'MtDyn$mtdyn'_'diff'_'$dyn'_'dyn'_'$kin'_'kin
#filename=R$cellRad'_'MtDyn$mtdyn'_'immo'_'$dyn'_'dyn'_'$kin'_'kin

## WIth cortical
#filename=R$cellRad'_'MtDyn$mtdyn'_'corticalDyn'_'$dyn'_'diff'_'$kin'_'kin
#filename=R$cellRad'_'MtDyn$mtdyn'_'diffEg5_surfDen_10_corticalDyn'_'$dyn'_'diff'_'$kin'_'kin'_'$lindensity
#filename=R$cellRad'_'MtDyn$mtdyn'_'corticalKin'_'$lindensity
#filename=R$cellRad'_'MtDyn$mtdyn'_'Eg5_surfDen10_corticalDynLinDen10'_'R$cellRad'_'DiffCoeff'_'$diffcoeff
#filename=R$cellRad'_'MtDyn$mtdyn'_'diffEg5_surfDen_10_corticalDyn_linDen_100_diffCoeff_20
filename=R$cellRad'_'MtDyn$mtdyn'_'Type'_'$mottype

echo $filename

#ffmpeg -r 60 -f image2 -s 1024x768 -i image%04d.ppm -vcodec libx264 -crf 25  -pix_fmt yuv420p  $filename'.avi'


#ffmpeg -r 60 -f image2 -i image%04d.ppm -vcodec libx264 -crf 25 $filename'.mp4'
#ffmpeg -f image2 -i image%04d.ppm -vcodec libx264 -crf 25 $filename'.mp4'
#ffmpeg -f image2  -r 10 -i image%04d.ppm -crf 25 $filename'.mp4'
#ffmpeg -i $filename'.mp4' -vcodec h264 -acodec mp2 $filename'2.mp4'


ffmpeg -f image2 -i image%04d.ppm -vcodec libx264 -crf 25 $filename'.mp4'
# ..... i need this
convert image0151.ppm $filename'_'final'.pdf'

#rm *.ppm
#rm result.out
#rm data.in
