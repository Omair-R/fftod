package tests

import "base:runtime"
import "core:math"
import "core:os"
import "core:c"

when os.OS == .Linux {
	foreign import fft_ooura "./vendor/fft4g.a"
} else when os.OS == .Windows {
	foreign import fft_ooura "./vendor/fft4g.lib"
}


@(default_calling_convention="c")
foreign fft_ooura {
	cdft :: proc(n: c.int, isgn: c.int, a: [^]f64, ip: [^]c.int, w: [^]f64) ---
}


Ooura_Plan :: struct
{
	n:int, 
	work: []f64,
	ip :[]c.int,
	twiddles: []f64,
}


make_ooura_plan :: proc(
	n: int,
) -> (
	plan: Ooura_Plan,
	err: runtime.Allocator_Error,
)#optional_allocator_error {

	sqrt_n := int(2+math.sqrt(f64(n)))

	work := make([]f64, 2*n) or_return
	twiddles := make([]f64, 2*n) or_return
	ip := make([]c.int, 2+sqrt_n)

	return Ooura_Plan{
		n,
		work,
		ip,
		twiddles,
	}, .None
}


delete_ooura_plan :: proc(
	plan: Ooura_Plan,
	allocator:=context.allocator,
) -> runtime.Allocator_Error {
	delete(plan.work)
	delete(plan.ip)
	delete(plan.twiddles)
	return .None
}


ooura_fft :: proc(x : []complex128, plan:Ooura_Plan) -> []complex128 {

	n := len(x)

	for i in 0..<n{
		plan.work[2*i] = real(x[i])
		plan.work[2*i+1] = imag(x[i])
	}
	cdft(
		2*c.int(n),
		1,
		raw_data(plan.work),
		raw_data(plan.ip),
		raw_data(plan.twiddles)
	)
	a := make([]complex128, n)

	for i in 0..<n{
		a[i] = complex(plan.work[2*i], plan.work[2*i+1])
	}

	return a
}

