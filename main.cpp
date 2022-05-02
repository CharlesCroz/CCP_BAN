#include <iostream>
#include "BigInt.h"

int main() {
    {
        unsigned char a_number[8] = {};
        a_number[0] = 0b00001111;
        a_number[3] = 0b00111100;
        a_number[5] = 0b11000011;
        a_number[7] = 0b00000011;

        BigInt<8> bi1(15);

        BigInt<8> bi2(a_number);
        std::cout << "bi1       " << bi1 << "\n";

        std::cout << "bi2       " << bi2 << "\n";
        BigInt<8> bi3(bi1 * bi2);
        std::cout << "bi3       " << bi3 << "\n";
        BigInt<8> bi4 = bi2 % bi1;
        std::cout << "bi2[bi1]  " << bi4 << "\n";
    }

    {
        BigInt<8> bi5(55555);
        BigInt<8> bi6(55555);
        BigInt<8> bi7(7);
        std::cout << "55555^55555[7]  " << BigInt<8>::pow_mod(bi5, bi6, bi7) << "\n";
    }

    {
        for (int i = 30; i < 40; i++) {
            try {
                std::cout << "2^" << i << " : " << BigInt<4>::power_of_two(i) << "\n";
            } catch (std::exception &e) {
                std::cout << e.what() << "\n";
                break;
            }
        }
    }

    {
        const char my_number[] = "99910801304194";
        BigInt<8> bi8(my_number);
        std::cout << my_number << "  " << bi8.to_string() << "\n";
    }

    {
        BigInt<8> bi9("123565498787564");
        bi9 %= BigInt<8>::power_of_two(18);
        std::cout << bi9 << "\n";
        std::cout << bi9.to_bytes() << "\n";

        BigInt<8>::byte_t dest[4];
        bi9.get_bytes(0, 4, dest);
        for (unsigned char i: dest) {
            std::cout << (int) i << " ";
        }
        std::cout << "\n";
        bi9.set_bytes(4, 8, dest);
        std::cout << bi9 << "\n";
        std::cout << bi9.to_bytes() << "\n";
    }

    {
        for (int i = 0; i < 5; ++i) {
            std::cout << BigInt64::rand() << "\n";
        }

        for (int i = 0; i < 5; ++i) {
            std::cout << BigInt64::rand(BigInt64::power_of_two(20 + i)) << "\n";
        }

        for (int i = 26; i < 30; ++i) {
            std::cout << BigInt64::rand(BigInt64::power_of_two(18 + i), BigInt64::power_of_two(20 + i)) << "\n";
        }
    }

    {
        BigInt32 bi32(0llu);
        std::cout << bi32.set_bit(0, true).set_bit(8, true).set_bit(16, true).set_bit(0, false) << "\n";

        std::cout << "4:" << bi32.get_bit(4) << ";8:" << bi32.get_bit(8) << "\n";
    }

    {
        BigInt64 bi_num("125478963");
        BigInt64 bi_denom("2458");
        std::cout << bi_num.to_bytes() << "\n";
        std::cout << bi_num << " / " << bi_denom << " = " << bi_num / bi_denom << "\n";
        std::cout << bi_num << " % " << bi_denom << " = " << bi_num % bi_denom << "\n";
        std::cout << BigInt64((bi_num / bi_denom) * bi_denom) + (bi_num % bi_denom) << "\n";
    }

    {
        BigInt128 bigPrime1("604541225532347786923774287563");
        std::cout << "Is " << bigPrime1 << " prime ? " << (bigPrime1.is_probably_prime(10) ? "Yes" : "No") << "\n\n";

        {
            BigInt256 bigPrime = BigInt256::generate_prime(256, 5);
            std::cout << "\n" << bigPrime << "\n is a big (" << bigPrime.get_size() << " bits) prime \n";
        }
        {
            BigInt512 bigPrime = BigInt512::generate_prime(512, 5);
            std::cout << "\n" << bigPrime << "\n is a big (" << bigPrime.get_size() << " bits) prime \n";
        }
    }

    return 0;
}
