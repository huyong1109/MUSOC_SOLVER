load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
filenames = (/"evp-solution", "bicg-solution","evp-bicg-sol"/)
do i = 0, dimsizes(filenames)-1
    filename = filenames(i) 
    path = filename+".bin"
    err = fbinrecread(path,0,(/92,182/),"double") 
    print(filename+" err min and max")
    printMinMax(err, True)
    
    wks   = gsn_open_wks ("png",filename);
    
    res                      = True                      ; plot mods desired
    res@cnFillOn             = True       ; turn on color fill
    res@cnLineLabelsOn       = False      ; turns off contour line labels
    res@cnFillMode           = "RasterFill"      ; turns off contour line labels
    ;res@cnLevelSelectionMode = "Explicitlevels"
    ;res@cnLevels = 10^ispan(-9,-2,1)
    ;res@cnLevels = fspan(-1,1,20)

    plot = gsn_csm_contour_map_ce(wks,err,res)     ; create plot
    system("convert -trim "+filename+".png " +filename+".png")
end do 

end
