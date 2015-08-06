load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
filename1 = "bicg-solution1"
filename = "bicg-solution"
path1 = filename1+".bin"
path = filename+".bin"
err1= fbinrecread(path1,0,(/92,182/),"double")
err = fbinrecread(path,0,(/92,182/),"double")
err2 = err1 - err
print("Err2 min and max")
printMinMax(err, True)

wks   = gsn_open_wks ("x11",filename);

res                   = True                      ; plot mods desired
  res@cnFillOn             = True       ; turn on color fill
  res@cnLineLabelsOn       = False      ; turns off contour line labels
  res@cnFillMode           = "RasterFill"      ; turns off contour line labels
  ;res@cnLevelSelectionMode = "Explicitlevels"
  ;res@cnLevels = 10^ispan(-9,-2,1)
  ;res@cnLevels = fspan(-1,1,20)

plot = gsn_csm_contour_map_ce(wks,err2,res)     ; create plot

end
