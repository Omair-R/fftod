package tests

import "core:slice"
import "core:fmt"
import "core:time"
import "core:math/rand"
import "core:strings"
import "core:testing"
import "core:log"
import fftod "../fftod"
import "core:math"


benchmark_setup :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error) {

	x := make([]complex128, options.bytes)
	for &val in x {
		val = complex(rand.float64(), rand.float64())
	}

	options.input = transmute([]u8)x
	return nil 
}


benchmark_teardown :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error) {

	x := transmute([]complex128)options.input
	delete(x)
	return nil
}


benchmark_print :: proc(
	str: ^strings.Builder,
	name: string,
	options: ^time.Benchmark_Options,
	loc := #caller_location
) {
	fmt.sbprintfln(str, "\n[%v - %v] processed \n%v values in %.4f ms %.1f r/s",
		name,
		options.bytes, 
		options.processed,
		time.duration_milliseconds(options.duration),
		options.rounds_per_second,
	)
}


Benchmark_Proc:: #type proc(
	_: ^time.Benchmark_Options,
	_:=context.allocator
) -> (err: time.Benchmark_Error)


benchmark_naive :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

	x := transmute([]complex128)options.input

	for _ in 0..=options.rounds{
		f := naive_dft(x)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}


benchmark_stockham_rad2 :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

		
	x := transmute([]complex128)options.input
	plan := fftod.make_stockham_plan(len(x), complex128)
	defer fftod.delete_stockham_plan(plan)

	for _ in 0..=options.rounds{
		f := fftod.out_fft_stockham_radix2(x, plan)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}


benchmark_stockham_rad4 :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

	x := transmute([]complex128)options.input
	plan := fftod.make_stockham_plan(len(x), complex128)
	defer fftod.delete_stockham_plan(plan)

	for _ in 0..=options.rounds{
		f := fftod.out_fft_stockham_radix4(x, plan)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}


benchmark_stockham_mixed :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

	x := transmute([]complex128)options.input
	plan := fftod.make_stockham_plan(len(x), complex128)
	defer fftod.delete_stockham_plan(plan)

	for _ in 0..=options.rounds{
		f := fftod.out_fft_stockham_mixed(x, plan)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}


benchmark_bluestein :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

	x := transmute([]complex128)options.input
	plan := fftod.make_bluestein_plan(len(x), complex128)
	defer fftod.delete_bluestein_plan(plan)

	for _ in 0..=options.rounds{
		f := fftod.out_fft_bluestein(x, plan)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}

benchmark_ooura :: proc(
	options: ^time.Benchmark_Options,
	allocator:=context.allocator
)->(err: time.Benchmark_Error){

	x := transmute([]complex128)options.input
	plan := make_ooura_plan(len(x))
	defer delete_ooura_plan(plan)

	for _ in 0..=options.rounds{
		f := ooura_fft(x, plan)
		delete(f)
	}
	options.count = options.rounds
	options.processed = options.rounds * options.bytes
	return nil
}


inner_benchmark :: proc(
	t: ^testing.T,
	str: ^strings.Builder,
	$rounds : int,
	signal_size : int,
) {

	bench_s := make([dynamic]string)
	bench_p := make([dynamic]Benchmark_Proc)

		
	// append(&bench_s, "Naive") // This is veeery slow
	// append(&bench_p, benchmark_naive)
	append(&bench_s, "Stockham Radix2")
	append(&bench_p, benchmark_stockham_rad2)
	append(&bench_s, "Stockham Radix4")
	append(&bench_p, benchmark_stockham_rad4)
	append(&bench_s, "Stockham Mixed")
	append(&bench_p, benchmark_stockham_mixed)
	append(&bench_s, "Bluestein")
	append(&bench_p, benchmark_bluestein)
	append(&bench_s, "Ooura")
	append(&bench_p, benchmark_ooura)

	for i in 0..<len(bench_s) {
		name := bench_s[i]
		options := &time.Benchmark_Options{
			rounds = rounds,
			bytes = signal_size,
			setup=benchmark_setup,
			teardown=benchmark_teardown,
			bench=bench_p[i],
		}
		err := time.benchmark(options)
		testing.expectf(t, err ==  nil, "Unexpected benchmark error : %v", name)
		if err == nil {
			benchmark_print(str, name, options)
		}		
	}

	delete(bench_s)
	delete(bench_p)
}


@(test)
benchmark_pow2 :: proc(t: ^testing.T) {
	str: strings.Builder
	strings.builder_init(&str, context.allocator)

	defer {
		log.info(strings.to_string(str))
		strings.builder_destroy(&str)
	}

	inner_benchmark(t, &str, 1_000, 64)
	inner_benchmark(t, &str, 1_000, 1024)
	// inner_benchmark(t, &str, 1_000, 4096)
}
