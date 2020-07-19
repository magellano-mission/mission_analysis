function const = PropagatedOrbits(data)

[T, Ysat] = ode113(@GaussIntegration, [0, data.tend], data.X0_kep, data.opt, data);

const.T = T;
const.Y = Ysat;

end