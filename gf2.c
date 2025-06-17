#include <stdexcept>
#include <iostream>

#include "gf2.h"

uint32_t gf2::bit_len(uint64_t num)
{
    uint32_t res = 0;

    while (num > 0) {
        num >>= 1;
        res++;
    }

    return res;
}

uint32_t gf2::reduce(uint64_t num, uint32_t mod)
{
    uint32_t num_len = bit_len(num);
    uint32_t mod_len = bit_len(mod);

    while (num_len >= mod_len)
    {
        num ^= (uint64_t)mod << (num_len - mod_len);
        num_len = bit_len(num);
    }

    return num;
}

bool gf2::operator==(const gf2 &other) const
{
    return _x == other._x && _mod == other._mod;
}

bool gf2::operator!=(const gf2 &other) const
{
    return _x != other._x || _mod != other._mod;
}

gf2 gf2::operator+ (const gf2 &other) const
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator+ _mod != other._mod");
    }

    gf2 res = other;
    res._x ^= _x;

    return res;
}

gf2 &gf2::operator+= (const gf2 &other)
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator+= _mod != other._mod");
    }

    _x ^= other._x;

    return *this;
}

gf2 gf2::operator- (const gf2 &other) const
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator- _mod != other._mod");
    }

    gf2 res = other;
    res._x ^= _x;

    return res;
}

gf2 &gf2::operator-= (const gf2 &other)
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator-= _mod != other._mod");
    }

    _x ^= other._x;

    return *this;
}

gf2 gf2::operator* (const gf2 &other) const
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator* _mod != other._mod");
    }

    uint64_t res = 0;
    uint32_t other_x = other._x;
    uint32_t deg = 0;

    while (other_x > 0)
    {
        res ^= (uint64_t)(_x * (other_x & 1)) << deg;
        other_x >>= 1;
        deg++;
    }

    return gf2(_mod, res);
}

gf2 &gf2::operator*= (const gf2 &other)
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator*= _mod != other._mod");
    }

    uint64_t res = 0;
    uint32_t other_x = other._x;
    uint32_t deg = 0;

    while (other_x > 0)
    {
        //std::cout << "other_x = " << other_x << std::endl;

        res ^= (uint64_t)(_x * (other_x & 1)) << deg;
        other_x >>= 1;
        deg++;
    }

    _x = reduce(res, _mod);

    return *this;
}

gf2 gf2::operator/ (const gf2 &other) const
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator/ _mod != other._mod");
    }

    uint32_t num = _x;
    uint32_t num_len = bit_len(num);
    uint32_t other_x_len = bit_len(other._x);
    uint32_t res = 0;

    while (num_len >= other_x_len)
    {
        res ^= 1 << (num_len - other_x_len);
        num ^= other._x << (num_len - other_x_len);
        num_len = bit_len(num);
    }

    return gf2(_mod, res);
}

gf2 &gf2::operator/= (const gf2 &other)
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator/= _mod != other._mod");
    }

    uint32_t num = _x;
    uint32_t num_len = bit_len(num);
    uint32_t other_x_len = bit_len(other._x);
    uint32_t res = 0;

    while (num_len >= other_x_len)
    {
        res ^= 1 << (num_len - other_x_len);
        num ^= other._x << (num_len - other_x_len);
        num_len = bit_len(num);
    }

    _x = res;

    return *this;
}

gf2 gf2::operator% (const gf2 &other) const
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator%% _mod != other._mod");
    }

    uint32_t res = reduce(_x, other._x);

    return gf2(_mod, res);
}

gf2 &gf2::operator%= (const gf2 &other)
{
    if (_mod != other._mod) {
        throw std::invalid_argument("gf2::operator%%= _mod != other._mod");
    }

    _x = reduce(_x, other._x);
    return *this;
}

gf2 gf2::operator^ (int32_t degree) const
{
    uint32_t mod = (1 << (bit_len(_mod)-1)) - 1;
    if (_x == 0 || _x == 1 || degree % mod == 1) {
        return *this;
    }
    if (degree % mod == 0) {
        return gf2(_mod, 1);
    }

    uint32_t pos_degree = degree % mod;
    if (degree < 0) {
        pos_degree = (int32_t)mod + (int32_t)(degree % (int32_t)mod);
    }

    gf2 res(_mod, 1);
    gf2 count = *this;

    while (pos_degree > 0) {
        res *= (pos_degree & 1) ? count : gf2(_mod, 1);
        pos_degree >>= 1;
        count *= count;
    }

    return res;
}

gf2 &gf2::operator^= (int32_t degree)
{
    *this = *this ^ degree;
 
    return *this;
}