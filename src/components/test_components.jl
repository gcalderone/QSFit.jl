
dom = Domain(1000.:1:8000)
iter=10_000

c = emline(4862)
test_component(c, dom, iter=iter)
c.profile = :Lorentzian
test_component(c, dom, iter=iter)

c = powerlaw(4862)
test_component(c, dom, iter=iter)

c = hostgalaxy("Ell5")
test_component(c, dom, iter=iter)

c = ironopt_broad(3000.)
test_component(c, dom, iter=iter)

c = ironopt_narrow(800.)
test_component(c, dom, iter=iter)

c = ironuv(30000.)
test_component(c, dom, iter=iter)


