fig = plt.figure()
ax = fig.add_subplot(111)
H, xedges, yedges = np.histogram2d(aa, bb, range=[[293.,1454.0], [464.,1896.0]], bins=(50, 50))
extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
subplots_adjust(bottom=0.15, left=0.15)
levels = (1.0e4, 1.0e3, 1.0e2, 2.0e1)
cset = contour(H, levels, origin=’lower’,colors=[‘black’,’green’,’blue’,’red’],linewidths=(1.9, 1.6, 1.5, 1.4),extent=extent)
plt.clabel(cset, inline=1, fontsize=10, fmt=’%1.0i’)
for c in cset.collections:
c.set_linestyle(‘solid’)