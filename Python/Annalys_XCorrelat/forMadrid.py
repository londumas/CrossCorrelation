
path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/XCorrelat_0_160_10Mpc_m2_Prod_DR12_11_CoAddDelta/source.root"
dr12_xi1D, dr12_xi2D, dr12_xi2D_2, dr12_xi2D_2_distX, dr12_xi2D_2_distY = get_Data(path,False)
path = "/home/gpfs/manip/mnt0607/bao/hdumasde/Results/RootFile/correlation_notDebug_allSky_eBOSS.root"
eBOSS_xi1D, eBOSS_xi2D, eBOSS_xi2D_2, eBOSS_xi2D_2_distX, eBOSS_xi2D_2_distY = get_Data(path,False)

fig = plt.figure()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.)

### Xi1D
coef = -1.
ax = fig.add_subplot(121)
plt.errorbar(dr12_xi1D[0,:],coef*dr12_xi1D[1,:],yerr=coef*dr12_xi1D[2,:],marker="o",label='DR12Q')
plt.errorbar(eBOSS_xi1D[0,:],coef*eBOSS_xi1D[1,:],yerr=coef*eBOSS_xi1D[2,:],marker="o",label='eBOSS')
plt.title(r'', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$\xi(|s|)$', fontsize=40)
myTools.deal_with_plot(False,False,True)

### Xi1D*|s|^{2}
coef = -1.*dr12_xi1D[0,:]*dr12_xi1D[0,:]
ax = fig.add_subplot(122)
plt.errorbar(dr12_xi1D[0,:],coef*dr12_xi1D[1,:],yerr=coef*dr12_xi1D[2,:],marker="o",label='DR12Q')
coef = -1.*eBOSS_xi1D[0,:]*eBOSS_xi1D[0,:]
plt.errorbar(eBOSS_xi1D[0,:],coef*eBOSS_xi1D[1,:],yerr=coef*eBOSS_xi1D[2,:],marker="o",label='eBOSS')
plt.title(r'', fontsize=40)
plt.xlabel(r'$|s| \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$|s|^{2}.\xi(|s|)$', fontsize=40)
myTools.deal_with_plot(False,False,True)
plt.show()


### Xi2D
coef = -1.
fig = plt.figure()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7, hspace=0.)
ax = fig.add_subplot(131)
ax.set_xticks([0.,50.,100.,150.])
plt.imshow(coef*dr12_xi2D_2, origin='lower',extent=[const.minY2D, const.maxY2D, const.minX2D, const.maxX2D],interpolation='None')
plt.title(r'$DR12Q: \, \xi(s_{\perp},s_{\parallel})$', fontsize=30, y=1.08)
plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
plt.ylabel(r'$s_{\parallel} \, [h^{-1} Mpc]$', fontsize=40)
plt.grid(True)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

### Xi2D*|s|^{2}
coef = -1.*numpy.sqrt(dr12_xi2D_2_distX*dr12_xi2D_2_distX + dr12_xi2D_2_distY*dr12_xi2D_2_distY)
ax = fig.add_subplot(132)
ax.set_xticks([0.,50.,100.,150.])
plt.imshow(coef*dr12_xi2D_2, origin='lower',extent=[const.minY2D, const.maxY2D, const.minX2D, const.maxX2D],interpolation='None') #interpolation='None'
plt.title(r'$DR12Q: \, |s|^{2}.\xi(s_{\perp},s_{\parallel})$', fontsize=30, y=1.08)
plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

### Xi2D*|s|^{2}
coef = -1.*numpy.sqrt(eBOSS_xi2D_2_distX*eBOSS_xi2D_2_distX + eBOSS_xi2D_2_distY*eBOSS_xi2D_2_distY)
ax = fig.add_subplot(133)
ax.set_xticks([0.,50.,100.,150.])
plt.imshow(coef*eBOSS_xi2D_2, origin='lower',extent=[const.minY2D, const.maxY2D, const.minX2D, const.maxX2D],interpolation='None') #interpolation='None'
plt.title(r'$eBOSS: \, |s|^{2}.\xi(s_{\perp},s_{\parallel})$', fontsize=30, y=1.08)
plt.xlabel(r'$s_{\perp} \, [h^{-1} Mpc]$', fontsize=40)
plt.grid(True)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

plt.show()


