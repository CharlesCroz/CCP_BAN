//
// Created by Charles on 25/02/2022.
//

#ifndef CPP_BAN_BIGINT_H
#define CPP_BAN_BIGINT_H

#include <string>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <stack>
#include <random>

#define BIGINT_BYTE_HEAD 0b10000000
#define BIGINT_BYTE_TAIL 0b00000001
#define BIGINT_BYTE_TWO  0b00000010
#define BIGINT_BYTE_ALL  0b11111111
#define BIGINT_BYTE_NONE 0b00000000
#define BIGINT_BYTE_MAX 0b100000000
#define BIGINT_BYTE_LEN 8
#define BIGINT_SEP "\'"
#define BIGINT_FCT to_string //or to_bytes

template<int BYTES>
class BigInt {
public:
    typedef unsigned char byte_t;
    typedef char16_t byte_long_t;
    typedef char32_t byte_longlong_t;
    typedef BigInt<BYTES> my_type;
    typedef BigInt<BYTES + 1> my_type_extended;
    typedef BigInt<BYTES * 2> my_type_doubled;
    typedef const my_type &const_ref;

    BigInt() = default;

    explicit BigInt(const byte_t val[BYTES]) {
        for (int i = 0; i < BYTES; ++i) {
            data[i] = val[i];
        }
    }

    template<int BYTES2>
    explicit BigInt(const BigInt<BYTES2> &other) {
        for (int i = 0; i < min(BYTES, BYTES2); ++i) {
            data[i] = other.get_byte(i);
        }
    }

    explicit BigInt(const char *words) {
        try {
            for (int i = 0; i < strlen(words); ++i) {
                if (words[i] < '0' || words[i] > '9') {
                    throw std::invalid_argument("argument must represent an integer");
                }
                this->times10();
                this->add_byte(words[i] - (byte_t) '0');
            }
        } catch (std::invalid_argument &e) {
            throw e;
        } catch (std::out_of_range &e) {
            throw std::invalid_argument("argument represents a number that is too big");
        }
    }

    BigInt(const BigInt &other) {
        for (int i = 0; i < BYTES; ++i) {
            data[i] = other.data[i];
        }
    }

    explicit BigInt(const unsigned long long int &val) {
        if (val > 0) {
            unsigned long long v = val;
            int i = 0;
            while (v > 0) {
                data[i] = v % BIGINT_BYTE_MAX;
                v = v / BIGINT_BYTE_MAX;
                ++i;
            }
        }
    }

    [[maybe_unused]]
    static my_type max_val() {
        my_type result;
        for (size_t i = 0; i < BYTES; ++i) {
            result.data[i] = BIGINT_BYTE_ALL;
        }
        return result;
    }

    static my_type rand() {
        my_type result = {};
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<unsigned char> distribution(0, BIGINT_BYTE_ALL);
        for (int i = 0; i < BYTES; ++i) {
            result.data[i] = distribution(gen);
        }
        return result;
    }

    static my_type rand(const_ref max) {
        my_type result = rand();
        result %= max;
        return result;
    }

    static my_type rand(const_ref min, const_ref max) {
        my_type result = rand();
        my_type diff = max - min;
        result %= diff;
        result += min;
        return result;
    }

    my_type &operator=(const_ref other) {
        if (this == &other) {
            return *this;
        }
        for (int i = 0; i < BYTES; ++i) {
            this->data[i] = other.data[i];
        }
        return *this;
    }

    [[maybe_unused]] [[nodiscard]]
    static my_type power_of_two(const int &i) {
        if (i >= (BYTES * BIGINT_BYTE_LEN)) {
            throw std::invalid_argument("Argument should not exceed 8*BYTES");
        } else if (i < 0) {
            throw std::invalid_argument("Argument should be a positive integer");
        }
        my_type result = {};
        result.data[i / BIGINT_BYTE_LEN] = 0b1 << (i % BIGINT_BYTE_LEN);
        return result;
    }

    [[maybe_unused]] [[nodiscard]]
    byte_t get_byte(const int &i) const {
        if (i > (BYTES * BIGINT_BYTE_LEN)) {
            throw std::invalid_argument("Argument must not exceed 8*BYTES");
        } else if (i < 0) {
            throw std::invalid_argument("Argument must be a positive integer");
        }
        return data[i];
    }

    // Returns bytes from i_min included to i_max excluded
    [[maybe_unused]]
    void get_bytes(const int &i_min, const int &i_max, byte_t *dest) const {
        if (i_min >= i_max) {
            throw std::invalid_argument("First argument must be lower than second argument");
        } else if (i_min < 0) {
            throw std::invalid_argument("First argument must be positive");
        } else if (i_max > BIGINT_BYTE_LEN * BYTES) {
            throw std::invalid_argument("Second argument must not exceed 8*BYTES");
        }
        for (int i = 0; i < (i_max - i_min); ++i) {
            dest[i] = data[i + i_min];
        }
    }

    [[maybe_unused]]
    my_type &set_byte(const int &i, const byte_t &b) {
        if (i < 0) {
            throw std::invalid_argument("first argument can't be negative");
        } else if (i > BIGINT_BYTE_LEN * BYTES - 1) {
            throw std::invalid_argument("first argument can't be bigger than the 8*BYTES");
        }
        data[i] = b;
        return *this;
    }

    my_type &set_bytes(const int &i_min, const int &i_max, const byte_t *src) {
        if (i_min >= i_max) {
            throw std::invalid_argument("First argument must be lower than second argument");
        } else if (i_min < 0) {
            throw std::invalid_argument("First argument must be positive or zero");
        } else if (i_max > BIGINT_BYTE_LEN * BYTES) {
            throw std::invalid_argument("Second argument must not exceed 8*BYTES");
        }
        for (int i = 0; i < (i_max - i_min); ++i) {
            data[i + i_min] = src[i];
        }
        return *this;
    }

    [[maybe_unused]] [[nodiscard]]
    bool get_bit(const int &i) const {
        return data[i / BIGINT_BYTE_LEN] & (BIGINT_BYTE_TAIL << (i % BIGINT_BYTE_LEN));
    }

    my_type &set_bit(const int &i, const bool &b) {
        int byte = i / BIGINT_BYTE_LEN;
        if (b) {
            data[byte] = data[byte] | (BIGINT_BYTE_TAIL << (i % BIGINT_BYTE_LEN));
        } else {
            data[byte] = data[byte] & ((BIGINT_BYTE_TAIL << (i % BIGINT_BYTE_LEN)) ^ BIGINT_BYTE_ALL);
        }
        return *this;
    }

    [[nodiscard]]
    std::string to_string() const {
        std::stack<int> digits;
        my_type me_too(*this);
        char decimal;
        while (me_too) {
            me_too = me_too.div10(decimal);
            digits.push(decimal);
        }
        std::string result = {};
        auto i = digits.size();
        if (i == 0) {
            result.append("0");
        } else {
            while (!digits.empty()) {
                result.append(1, (char) ('0' + digits.top()));
                digits.pop();
                --i;
                if (i % 3 == 0 && i > 0) {
                    result.append(BIGINT_SEP);
                }
            }
        }
        return result;
    }

    [[maybe_unused]] [[nodiscard]]
    std::string to_bytes() const {
        std::string result = "[";
        char buffer[10];
        for (int i = 0; i < BYTES; ++i) {
            sprintf(buffer, "%3d", data[i]);
            result.append(buffer);
            if (i != BYTES - 1) {
                result.append(" ");
            } else {
                result.append("]");
            }
        }
        return result;
    }

    my_type &times2() {
        bool carry = false;
        int tmp;
        for (int i = 0; i < BYTES; ++i) {
            tmp = data[i];
            data[i] = tmp << 1;
            if (carry) {
                data[i] = data[i] | BIGINT_BYTE_TAIL;
            }
            carry = (tmp & BIGINT_BYTE_HEAD) == BIGINT_BYTE_HEAD;
        }
        return *this;
    }

    my_type &div2() {
        bool carry = false;
        int tmp;
        for (int i = BYTES - 1; i >= 0; --i) {
            tmp = data[i];
            data[i] = tmp >> 1;
            if (carry) {
                data[i] = data[i] | BIGINT_BYTE_HEAD;
            }
            carry = (tmp & BIGINT_BYTE_TAIL) == BIGINT_BYTE_TAIL;
        }
        return *this;
    }

    my_type &times10() {
        byte_long_t tmp;
        byte_long_t carry = 0;
        for (int i = 0; i < BYTES; ++i) {
            tmp = this->data[i] * 10 + carry;
            data[i] = (byte_t) (tmp & BIGINT_BYTE_ALL);
            carry = tmp >> BIGINT_BYTE_LEN;
        }
        if (carry > 0) {
            throw std::out_of_range("times10 overflow");
        }
        return *this;
    }

    my_type &div10(char &res) {
        byte_long_t carry = 0;
        byte_long_t tmp;
        for (int i = BYTES - 1; i >= 0; --i) {
            tmp = (data[i] + (carry << BIGINT_BYTE_LEN));
            data[i] = (byte_t) (tmp / 10);
            carry = tmp % 10;
        }
        res = (char) carry;
        return *this;
    }

    my_type &add_byte(const byte_t &b) {
        byte_long_t tmp = b;
        int i = 0;
        while (tmp != 0 && i < BYTES * BIGINT_BYTE_LEN) {
            tmp += data[i];
            data[i] = tmp & BIGINT_BYTE_ALL;
            tmp = tmp >> BIGINT_BYTE_LEN;
            ++i;
        }
        return *this;
    }

    my_type &mult_mod(const_ref other, const_ref m) {
        my_type_doubled product = *this * other;
        product %= my_type_doubled(m);
        for (int i = 0; i < BYTES; ++i) {
            data[i] = product.get_byte(i);
        }
        return *this;
    }

    my_type &mult_mod(const_ref other, const my_type_doubled &m) {
        my_type_doubled product = *this * other;
        product %= m;
        for (int i = 0; i < BYTES; ++i) {
            data[i] = product.get_byte(i);
        }
        return *this;
    }

    [[maybe_unused]] [[nodiscard]]
    static my_type pow_mod(const_ref base,
                           const_ref exp,
                           const_ref m) {
        if (m.is_zero()) {
            throw std::invalid_argument("third argument can't be zero");
        }
        my_type result = my_type::one();
        my_type exp_too(exp);
        my_type_doubled big_m(m);
        my_type growing(base % m);
        while (exp_too.not_zero()) {
            if (exp_too.is_odd()) {
                result.mult_mod(growing, big_m);
            }
            growing.mult_mod(growing, big_m);
            exp_too.div2();
        }

        return result;
    }

    [[maybe_unused]] [[nodiscard]]
    static my_type one() {
        byte_t data[BYTES] = {};
        data[0] = 1;
        return my_type(data);
    }

    my_type &operator|=(const_ref other) {
        for (int i = 0; i < BYTES; ++i) {
            data[i] = data[i] | other.data[i];
        }
    }

    my_type operator|(const_ref other) {
        return my_type(*this).operator!=(other);
    }

    my_type &operator&=(const_ref other) {
        for (int i = 0; i < BYTES; ++i) {
            data[i] = data[i] & other.data[i];
        }
        return *this;
    }

    my_type operator&(const_ref other) {
        return my_type(*this).operator&=(other);
    }

    my_type &operator^=(const_ref other) {
        for (int i = 0; i < BYTES; ++i) {
            data[i] = data[i] ^ other.data[i];
        }
        return *this;
    }

    my_type operator^(const_ref other) {
        return my_type(*this).operator^=(other);
    }

    my_type &operator++() {
        int i = 0;
        bool carry = true;
        while (i < BYTES && carry) {
            if (data[i] < BIGINT_BYTE_ALL) {
                ++data[i];
                carry = false;
            } else {
                data[i] = 0;
            }
            ++i;
        }
        return *this;
    }

    my_type operator++(int) {
        my_type that(*this);
        ++*this;
        return that;
    }

    my_type &operator--() {
        int i = 0;
        bool carry = true;
        while (i < BYTES && carry) {
            if (data[i] > 0) {
                --data[i];
                carry = false;
            } else {
                data[i] = BIGINT_BYTE_ALL;
            }
            ++i;
        }
        return *this;
    }

    my_type operator--(int) {
        my_type that(*this);
        --*this;
        return that;
    }

    [[nodiscard]]
    my_type operator-(const_ref other) const {
        my_type result;
        int carry = 0;
        int tmp;
        for (int i = 0; i < BYTES; ++i) {
            tmp = this->data[i] - other.data[i] - carry;
            if (tmp < 0) {
                tmp += BIGINT_BYTE_MAX;
                carry = 1;
            } else {
                carry = 0;
            }
            result.data[i] = (byte_t) tmp;
        }
        return result;
    }

    void operator-=(const_ref other) {
        int carry = 0;
        int tmp;
        for (int i = 0; i < BYTES; ++i) {
            tmp = this->data[i] - other.data[i] - carry;
            if (tmp < 0) {
                tmp += BIGINT_BYTE_MAX;
                carry = 1;
            } else {
                carry = 0;
            }
            this->data[i] = (byte_t) tmp;
        }
    }

    my_type &operator%=(const_ref other) {
        if (other.is_zero()) {
            throw std::invalid_argument("Can't divide by zero");
        } else if (*this == other) {
            *this = {};
            return *this;
        }

        int shifts = 0;
        my_type_extended me_too(*this);
        my_type_extended tmp(other);
        my_type_extended next_tmp(other);
        next_tmp.times2();

        while (me_too > next_tmp) {
            tmp = next_tmp;
            next_tmp = next_tmp.times2();
            ++shifts;
        }
        for (int i = 0; i <= shifts; ++i) {
            if (me_too >= tmp) {
                me_too -= (tmp);
            }
            tmp.div2();
        }
        for (int i = 0; i < BYTES; ++i) {
            data[i] = me_too.get_byte(i);
        }
        return *this;
    }

    my_type operator%(const_ref other) const {
        return my_type(*this).operator%=(other);
    }

    my_type operator/(const_ref other) {
        return my_type(*this).operator/=(other);
    }

    [[nodiscard]]
    my_type &operator/=(const_ref other) {
        if (other.is_zero()) {
            throw std::invalid_argument("Can't divide by zero");
        }

        int shifts = 0;
        my_type_extended me_too(*this);
        my_type_extended tmp(other);
        my_type_extended next_tmp(other);
        next_tmp.times2();

        *this = {};

        while (me_too > next_tmp) {
            tmp = next_tmp;
            next_tmp = next_tmp.times2();
            ++shifts;
        }
        for (int i = 0; i <= shifts; ++i) {
            if (me_too >= tmp) {
                me_too -= tmp;
                this->set_bit(shifts - i, true);
            }
            tmp.div2();
        }
        return *this;
    }

    [[nodiscard]]
    my_type_doubled operator*(const_ref other) const {
        byte_longlong_t long_res[2 * BYTES] = {};
        byte_t res_data[2 * BYTES];
        for (int i = 0; i < BYTES; ++i) {
            for (int j = 0; j < BYTES; j++) {
                long_res[i + j] += this->data[i] * other.data[j];
            }
        }
        byte_longlong_t carry = 0;
        for (int i = 0; i < 2 * BYTES; ++i) {
            long_res[i] += carry;
            res_data[i] = (byte_t) (long_res[i] & BIGINT_BYTE_ALL);
            carry = long_res[i] >> BIGINT_BYTE_LEN;
        }
        return my_type_doubled(res_data);
    }

    [[nodiscard]]
    my_type operator+(const_ref other) const {
        my_type result;
        int carry = 0;
        int tmp;
        for (int i = 0; i < BYTES; ++i) {
            tmp = carry + (int) this->data[i] + (int) other.data[i];
            carry = tmp >> 8;
            result.data[i] = tmp & BIGINT_BYTE_ALL;
        }
        return result;
    }

    void operator+=(const_ref other) {
        int carry = 0;
        int tmp;
        for (int i = 0; i < BYTES; ++i) {
            tmp = carry + (int) this->data[i] + (int) other.data[i];
            carry = tmp >> 8;
            data[i] = tmp & BIGINT_BYTE_ALL;
        }
    }

    bool operator==(const_ref other) const {
        if (&other == this) {
            return true;
        }
        for (int i = 0; i < BYTES; ++i) {
            if (this->data[i] != other.data[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const_ref other) const {
        return !this->operator==(other);
    }

    bool operator>(const_ref other) const {
        if (&other == this) {
            return false;
        }
        for (int i = BYTES - 1; i >= 0; --i) {
            if (this->data[i] > other.data[i]) {
                return true;
            } else if (this->data[i] < other.data[i]) {
                return false;
            }
        }
        return false;
    }

    bool operator>=(const_ref other) const {
        if (&other == this) {
            return true;
        }
        for (int i = BYTES - 1; i >= 0; --i) {
            if (this->data[i] > other.data[i]) {
                return true;
            } else if (this->data[i] < other.data[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator<(const_ref other) const {
        return !this->operator>=(other);
    }

    bool operator<=(const_ref other) const {
        return !this->operator>(other);
    }

    explicit operator bool() const {
        for (int i = 0; i < BYTES; ++i) {
            if (data[i] != 0) {
                return true;
            }
        }
        return false;
    }

    [[maybe_unused]] [[nodiscard]]
    bool is_zero() const {
        for (int i = 0; i < BYTES; ++i) {
            if (data[i] != 0) {
                return false;
            }
        }
        return true;
    }

    [[maybe_unused]] [[nodiscard]]
    bool is_one() const {
        if (data[0] != BIGINT_BYTE_TAIL) {
            return false;
        }
        for (int i = 1; i < BYTES; ++i) {
            if (data[i] != 0) {
                return false;
            }
        }
        return true;
    }

    [[maybe_unused]] [[nodiscard]]
    bool is_two() const {
        if (data[0] != BIGINT_BYTE_TWO) {
            return false;
        }
        for (int i = 1; i < BYTES; ++i) {
            if (data[i] != 0) {
                return false;
            }
        }
        return true;
    }

    [[maybe_unused]] [[nodiscard]]
    bool not_zero() const {
        return !is_zero();
    }

    [[maybe_unused]] [[nodiscard]]
    bool is_odd() const {
        return (data[0] & BIGINT_BYTE_TAIL) == 1;
    }

    [[maybe_unused]] [[nodiscard]]
    bool is_even() const {
        return (data[0] & BIGINT_BYTE_TAIL) == 0;
    }

    [[maybe_unused]][[nodiscard]]
    bool is_probably_prime(const size_t tests) {
        if (is_zero() || is_one()) {
            return false;
        }
        if (is_two()) {
            return true;
        } else if (is_even()) {
            return false;
        }

        my_type n_minus_one(*this);
        --n_minus_one;

        my_type d(n_minus_one);
        size_t s = 0;

        while (d.is_even()) {
            ++s;
            d.div2();
        }

        bool is_composite;
        for (size_t i = 0; i < tests; ++i) {
            is_composite = true;
            my_type x = my_type::pow_mod(my_type::rand(my_type::power_of_two(1), n_minus_one), d, *this);
            if (x.is_one() || (x == n_minus_one)) {
                is_composite = false;
            } else {
                for (size_t j = 0; j < s; ++j) {
                    x.mult_mod(x, *this);
                    if (x == n_minus_one) {
                        is_composite = false;
                        break;
                    }
                }
            }

            if (is_composite) {
                return false;
            }
        }

        return true;
    }

    [[maybe_unused]][[nodiscard]]
    static my_type generate_prime_candidate(const unsigned int bits) {
        if (bits < 1) {
            throw std::invalid_argument("Argument must be strictly positive");
        } else if (bits > BIGINT_BYTE_LEN * BYTES) {
            throw std::invalid_argument("Argument must be lower than 8 * BYTES");
        }
        my_type candidate = my_type::rand();
        candidate.set_bit(0, true);
        candidate.set_bit(bits - 1, true);
        for (size_t i = bits; i < BIGINT_BYTE_LEN * (bits / BIGINT_BYTE_LEN); ++i) {
            candidate.set_bit(i, false);
        }
        for (size_t i = (bits / BIGINT_BYTE_LEN) + 1; i < BYTES; ++i) {
            candidate.data[i] = BIGINT_BYTE_NONE;
        }
        return candidate;
    }

    [[maybe_unused]][[nodiscard]]
    static my_type generate_prime(const unsigned int bits, const size_t tests) {
        my_type candidate;
        do {
            candidate = my_type::generate_prime_candidate(bits);
        } while (!candidate.is_probably_prime(tests));
        return candidate;
    }

    [[maybe_unused]][[nodiscard]]
    size_t get_size() const {
        return BIGINT_BYTE_LEN * BYTES;
    }

private:

    byte_t data[BYTES] = {};

    int min(int a, int b) {
        return (a < b ? a : b);
    }

};

[[maybe_unused]] typedef BigInt<4> BigInt32;
[[maybe_unused]] typedef BigInt<8> BigInt64;
[[maybe_unused]] typedef BigInt<16> BigInt128;
[[maybe_unused]] typedef BigInt<32> BigInt256;
[[maybe_unused]] typedef BigInt<64> BigInt512;
[[maybe_unused]] typedef BigInt<128> BigInt1024;
[[maybe_unused]] typedef BigInt<256> BigInt2048;

template<int BYTES>
[[maybe_unused]][[nodiscard]]
BigInt<BYTES> min(const BigInt<BYTES> &a, const BigInt<BYTES> &b) {
    return (a > b ? b : a);
}

template<int BYTES>
[[maybe_unused]][[nodiscard]]
BigInt<BYTES> max(const BigInt<BYTES> &a, const BigInt<BYTES> &b) {
    return (a > b ? a : b);
}

template<int BYTES>
std::ostream &operator<<(std::ostream &os, const BigInt<BYTES> &val) {
    os << val.BIGINT_FCT();
    return os;
}

#endif //CPP_BAN_BIGINT_H
