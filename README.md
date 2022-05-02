# BigAssNumber

Header Only class to handle (most likely quite poorly) cryptographic numbers.

## Requirements

C++17

## How To

Copy the .h file somewhere. Include it. Profit. 

Class is templated over how many bytes you want.

`BigInt<8> myBigInt(42)` for 42 as an 8 bytes (64 bits) integer.

`BigInt64 myBigInt(42)` is an alias.

`BigInt<16> myBigInt("12345678901234567891234567890123456789")`

`BigInt<32> myBigInt=BigInt<32>::power_of_two(200)`

## Nice things

`BigInt<BYTES>::pow_mod(a, n, m)` to compute `a^n[m]`
with a, n and m of type `BigInt<BYTES>`.

`a.mult_mod(b, m)` to compute `(a * b) % m`

`BigInt64::rand()`, `BigInt64::rand(max)`, `BigInt64::rand(min, max)`

## Untemplated classes

`BigIng64`, `BigIng128`, `BigIng256`, `BigIng512` as aliases for
`BigInt<8>`, `BigInt<16>`, `BigInt<32>`, `BigInt<64>` 



