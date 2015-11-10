package main

import (
	"fmt"
	//"log"
	"math"
	"os"
)

// coefficient and general approach adopted from:
//Adapted from R. Sonntag, C. Borgnakke and G. J. Wylen,
// Fundamentals of Classical Thermodynamics,
//5th Ed. John Wiley & Sons, 1998.
const (
	// universal gas constant:
	R        float64 = 8.314
	wr       float64 = 0.3978
	atmtobar float64 = 1.01325
	// s denotes simple molecules (methane)
	a1s    float64 = 0.1181193
	a2s    float64 = 0.265728
	a3s    float64 = 0.154790
	a4s    float64 = 0.030323
	b1s    float64 = 0.0236744
	b2s    float64 = 0.0186984
	b3s    float64 = 0
	b4s    float64 = 0.042724
	c1s    float64 = 0.155488e-4
	c2s    float64 = 0.623689e-4
	betas  float64 = 0.65392
	gammas float64 = 0.060167
	// r denotes reference molecules (n-octane)
	a1r    float64 = 0.2026579
	a2r    float64 = 0.331511
	a3r    float64 = 0.027655
	a4r    float64 = 0.203488
	b1r    float64 = 0.0313885
	b2r    float64 = 0.0503618
	b3r    float64 = 0.016901
	b4r    float64 = 0.041577
	c1r    float64 = 0.48736e-4
	c2r    float64 = 0.740336e-5
	betar  float64 = 1.226
	gammar float64 = 0.03754
)

type Result struct {
	Z, Z0, Z1    float64 // compressibility factors (total: Z, simple: 0...)
	Hr, Hr0, Hr1 float64 // residual enthalpies (total: Hr, simple: 0...)
	Sr, Sr0, Sr1 float64 // residual entropies (total: Sr, simple: 0...)
	F, F0, F1    float64 // fugacity coefficient
}

func (res Result) Printout() {
	fmt.Println("This is the table!")
	fmt.Println("\t Total \tSimple \t Correction")
	fmt.Printf("Z: \t%.4f \t%.4f\t %.4f\n", res.Z, res.Z0, res.Z1)
	fmt.Printf("f: \t%.4f \t%.4f\t %.4f\n", res.F, res.F0, res.F1)
	fmt.Printf("Hr: \t%.4f \t%.4f\t %.4f\n", res.Hr, res.Hr0, res.Hr1)
	fmt.Printf("Sr: \t%.4f \t%.4f\t %.4f\n", res.Sr, res.Sr0, res.Sr1)
}

func main() {

	var (
		Tc    float64 = 370 // Kelvin
		Pc    float64 = 42  // bar
		Tboil float64 = 230 // Kelvin
		T     float64 = 300 // Kelvin
		P     float64 = 10  // bar
		//w     float64 = 0.150
	)
	w := acentric_factor(Tboil, Tc, Pc)
	res := LeeKesler(Tc, Pc, T, P, w)
	res.Printout()
}
func LeeKesler(Tc, Pc, T, P, w float64) (r Result) {
	_ = w
	var res Result // result variable to be filled and returned
	Tr := T / Tc
	Pr := P / Pc
	// We'll use Vguess from ideal gas law (or uncomment and overwrite!)
	//Vguess := (R * T) / (P * 100000)
	var Vguess float64 = 4
	fmt.Println("Reduced State: Tr, Pr, Vguess: ", Tr, Pr, Vguess)
	// simple:
	As := a1s - a2s/Tr - a3s/(Tr*Tr) - a4s/math.Pow(Tr, 3)
	Bs := b1s - b2s/Tr + b3s/math.Pow(Tr, 3)
	Cs := c1s + c2s/Tr
	// units??
	EOSs := func(Vr float64) float64 {
		return (Tr / Vr) * (1 + (As / Vr) + Bs/(Vr*Vr) + Cs/math.Pow(Vr, 5) +
			(b4s/(math.Pow(Tr, 3)*Vr*Vr))*(betas+gammas/math.Pow(Vr, 2))*math.Exp(-gammas/math.Pow(Vr, 2)))
	}
	simple := func(Vr float64) float64 {
		return EOSs(Vr) - Pr
	}
	//////////////////////////////////////////////////////////////////////
	// reference:
	Ar := a1r - a2r/Tr - a3r/(Tr*Tr) - a4r/math.Pow(Tr, 3)
	Br := b1r - b2r/Tr + b3r/math.Pow(Tr, 3)
	Cr := c1r + c2r/Tr
	// units??
	EOSr := func(Vr float64) float64 {
		return (Tr / Vr) * (1 + (Ar / Vr) + Br/(Vr*Vr) + Cr/math.Pow(Vr, 5) +
			(b4r/(math.Pow(Tr, 3)*Vr*Vr))*(betar+gammar/math.Pow(Vr, 2))*math.Exp(-gammar/math.Pow(Vr, 2)))
	}
	reference := func(Vr float64) float64 {
		return EOSr(Vr) - Pr
	}
	//PlotToFile(EOSr, 0.00001, 15.0000, 0.001, "ref2")
	//PlotToFile(simple, 0.0001, 15.0000, 0.001, "ref3")
	Vrs := Rootfinder(simple, Vguess, 0.00001)
	Vrr := Rootfinder(reference, Vguess, 0.00001)
	fmt.Println("Vrs = ", Vrs)
	fmt.Println("Vrr = ", Vrr)
	Zref := (Pr * Vrr) / Tr
	//fmt.Println("Zref je:", Zref)
	res.Z0 = (Pr * Vrs) / Tr
	res.Z1 = (Zref - res.Z0) / wr
	res.Z = res.Z0 + w*res.Z1
	// that was compressibility, let's try fugacity now:
	Ds := b4s / (2 * math.Pow(Tr, 3) * gammas) *
		(betas + 1 -
			(betas+1+gammas/math.Pow(Vrs, 2))*
				math.Exp(-gammas/math.Pow(Vrs, 2)))

	lnFugStd := res.Z0 - 1 - math.Log(res.Z0) + As/Vrs +
		Bs/(2*Vrs*Vrs) +
		Cs/(5*math.Pow(Vrs, 5)) + Ds
	///////////////////////////////////////////////////
	Dr := b4r / (2 * math.Pow(Tr, 3) * gammar) *
		(betar + 1 -
			(betar+1+gammar/math.Pow(Vrr, 2))*
				math.Exp(-gammar/math.Pow(Vrr, 2)))
	lnFugRef := Zref - 1 - math.Log(Zref) + Ar/Vrr +
		Br/(2*Vrr*Vrr) +
		Cr/(5*math.Pow(Vrr, 5)) + Dr
	res.F0 = math.Exp(lnFugStd)
	res.F1 = math.Exp((lnFugRef - lnFugStd) / wr)
	lntot := lnFugStd + (w/wr)*(lnFugRef-lnFugStd)
	fmt.Println(lntot)
	//////////////////////////////////////////////////
	// Next: enthalpy departure:
	Hs_dep := Tr *
		(res.Z0 - 1 -
			(a2s+(2*a3s/Tr)+3*a4s/math.Pow(Tr, 2))/(Tr*Vrs) -
			(b2s-3*b3s/math.Pow(Tr, 2))/(2*Tr*math.Pow(Vrs, 2)) +
			c2s/(5*Tr*math.Pow(Vrs, 5)) +
			3*Ds)
	Hr_dep := Tr *
		(Zref - 1 -
			(a2r+(2*a3r/Tr)+3*a4r/math.Pow(Tr, 2))/(Tr*Vrr) -
			(b2r-3*b3r/math.Pow(Tr, 2))/(2*Tr*math.Pow(Vrr, 2)) +
			c2r/(5*Tr*math.Pow(Vrr, 5)) +
			3*Dr)
	//fmt.Println(Hs_dep, Hr_dep)
	res.Hr0 = Hs_dep
	res.Hr1 = (Hr_dep - res.Hr0) / wr
	res.Hr = res.Hr0 + w*res.Hr1
	//fmt.Println("\tA\t\t\t B\t\t\t C\t\tD")
	//fmt.Println(As, Bs, Cs, Ds)
	//fmt.Println(Ar, Br, Cr, Dr)
	/// Finally enthropy departure:
	//  ???????????????? where is this gone??????math.Log(P/atmtobar) -
	Ss_dep := math.Log(res.Z0) -
		(a1s+a3s/math.Pow(Tr, 2)+2*a4s/math.Pow(Tr, 3))/Vrs -
		(b1s-2*b3s/math.Pow(Tr, 3))/(2*math.Pow(Vrs, 2)) -
		c1s/(5*math.Pow(Vrs, 5)) +
		2*Ds
	Sr_dep := math.Log(Zref) -
		(a1r+a3r/math.Pow(Tr, 2)+2*a4r/math.Pow(Tr, 3))/Vrr -
		(b1r-2*b3r/math.Pow(Tr, 3))/(2*math.Pow(Vrr, 2)) -
		c1r/(5*math.Pow(Vrr, 5)) +
		2*Dr
	res.Sr0 = Ss_dep
	res.Sr1 = (Sr_dep - res.Sr0) / wr
	res.Sr = res.Sr0 + w*res.Sr1

	return res
}

func acentric_factor(Tb, Tc, Pc float64) float64 {
	// assume boiling pressure of 1 atm == 1.01325 bar (+ units are in bar and Kelvin)
	Pbr := atmtobar / Pc
	Tbr := Tb / Tc
	w := (math.Log(Pbr) - 5.92714 + 6.09648/Tbr + 1.28862*math.Log(Tbr) - 0.169347*math.Pow(Tbr, 6)) /
		(15.2518 - 15.6875/Tbr - 13.4721*math.Log(Tbr) + 0.43577*math.Pow(Tbr, 6))
	return w
}
func Rootfinder(f func(X float64) float64, guess, tolerance float64) float64 {
	// fine tune this if getting weird results:
	var delta float64 = 0.0001
	var maxiter int = 100000
	//pseudoderivative:
	df := func(x float64) float64 {
		dx := delta
		return ((f(x+dx) - f(x)) / dx)
	}
	//initialization
	x := guess
	x = x - f(guess)/df(guess)
	// main loop, using pseudoderivative
	i := 0
	for math.Abs(f(x)) > tolerance {
		x = x - f(x)/df(x)
		i++
		if i > maxiter {
			break
		}
	}
	//fmt.Println("error of measurement: ", f(x))
	fmt.Println("Derivative at f(x) ~ 0:", df(x))
	return x

}
func PlotToFile(f func(X float64) float64, start, end, step float64, filename string) {
	var vsetko string
	for i := start; i <= end; i = i + step {
		vsetko = vsetko + fmt.Sprintf("%g %g\n", i, f(i))

	}
	file, err := os.Create(filename)
	if err != nil {
		//	log.Fatal(err)
	}
	defer file.Close()
	file.WriteString(vsetko)
}
