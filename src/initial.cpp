// initial data generation that is incomplete. 
// for now (3/2/2022) we just use Mathematica to generate initial stuff, but ideally this should be done here instead so we may more easily experiment with setups.
// im just moving it here becuase i dont want it polluting numerical.cpp


//
//Initial Condition Functions
//

f64 P(f64 t) {
	if (t < t0 - delta / 2) return -k / 2;
	if (t > t0 + delta / 2) return  k / 2;
	f64 inside = (t - t0) / (sqrt((t - t0 + delta / 2) * (t0 - t + delta / 2)));
	return k / 2 * tanh(inside);
}

array<f64> PSample(f64 dt) {
	array<f64> values;
	for (int i = 0; i < Nt; i++) values.add(P(i * dt));
	return values;
}

f64 initialPos(f64 rho) {
	f64 M0 = mmm / (3 * M_SQRT_TWO * gg);
	f64 gamma = sqrt(1 + (k / 2) * (k / 2) / (M0 * M0));
	return m / (2 * g) * tanh(gamma * m * rho / M_SQRT_TWO);
}

array<f64> initialPosSample() {
	f64 dx = L / (2 * Nt + 1);
	array<f64> values;
	for (int i = 0; i < 2 * Nt + 1; i++) values.add(initialPos(-L / 2 + i * dx));
	return values;
}

f64 IIRiemannSum(f64 hk, f64 M0, f64 gamma) {
	f64 sum = 0;
	f64 dx = 1 / (100 * m);
	for (f64 i = -L * 100; i < L * 100; i += dx) {
		f64 s0 = sech(i);
		f64 s1 = sech(gamma * i);
		sum += s0 * s0 * s1 * s1 * dx;
	}
	return sum;
}

f64 initialMomentum(f64 rho) {
	f64 hk = k / 2;
	f64 M0 = mmm / (3 * M_SQRT_TWO * gg);
	f64 gamma = sqrt(1 + hk * hk / (M0 * M0));
	f64 II = IIRiemannSum(hk, M0, gamma);
	f64 sech0 = sech(sqrt(1 + hk * hk / (M0 * M0) * m * rho / M_SQRT_TWO));
	f64 sech1 = sech(m * rho / M_SQRT_TWO);
	return 3 * g * k / (4 * m) * ( sech0 * sech0 - (3 * II * sech1 * sech1) / 4);
}

array<f64> initialMomentumSample() {
	f64 dx = L / (2 * Nt + 1);
	array<f64> values;
	for (int i = 0; i < 2 * Nt + 1; i++) values.add(initialMomentum(-L / 2 + i * dx));
	return values;
}

array<f64> psinotSample() {
	f64 dx = L / (2 * Nt + 1);
	f64 coeff = sqrt(3 * m / (4 * (M_SQRT_TWO)));
	array<f64> values;
	for (int i = 0; i < 2 * Nt; i++) {
		f64 sh = sech(m / M_SQRT_TWO * (-L / 2 + i * dx));
		values.add(coeff * sh * sh);
	}
	return values;
}

array<f64> dpsinotSample() {
	f64 dx = L / (2 * Nt + 1);
	f64 coeff = -sqrt(3 * mmm / (2 * (M_SQRT_TWO)));
	array<f64> values;
	for (int i = 0; i < 2 * Nt; i++) {
		f64 inner = m / M_SQRT_TWO * (-L / 2 + i * dx);
		f64 sh = sech(inner);
		f64 th = tanh(inner);
		values.add(coeff * sh * sh * th);
	}
	return values;
}
