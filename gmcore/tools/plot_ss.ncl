f = addfile("../run/GMCORE-TESTBED/ss.360x180/ss.360x180x26.h0.nc", "r")

wks = gsn_open_wks("pdf", "ss.360x180x26")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@cnLinesOn = False
res@cnFillOn = True
res@trYReverse = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 35
res@cnLevelSpacingF = 5

u_zavg = zonalAve(f->u)

res@tiMainString = "0 Day"
plot_0  = gsn_csm_contour(wks, u_zavg(0,:,:), res)
res@tiMainString = "30 Day"
plot_30 = gsn_csm_contour(wks, u_zavg(30,:,:), res)

res := True

gsn_panel(wks, (/plot_0,plot_30/), (/1,2/), res)
