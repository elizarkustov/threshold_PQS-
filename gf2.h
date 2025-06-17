#ifndef GF2_H
#define GF2_H

#include <cstdint>

class gf2
{
private:
    uint32_t _x;
    uint32_t _mod;
public:
    static uint32_t bit_len(uint64_t num);
    static uint32_t reduce(uint64_t num, uint32_t mod);

    gf2(): _x(0), _mod(1) {}
    gf2(uint32_t mod): _x(0), _mod(mod) {}
    gf2(uint32_t mod, uint64_t x): _x(0), _mod(mod) {
        _x = reduce(x, _mod);
    }
    gf2(const gf2& other): _x(other._x), _mod(other._mod) {}

    gf2 &operator=(uint64_t x)
    {
        _x = reduce(x, _mod);

        return *this;
    }
    gf2 &operator=(const gf2 &other)
    {
        _mod = other._mod;
        _x = other._x;

        return *this;
    }

    operator uint32_t() const {return _x;}
    
    bool operator==(const gf2 &other) const;
    bool operator!=(const gf2 &other) const;

    gf2 operator+ (const gf2 &other) const;
    gf2 &operator+= (const gf2 &other);
    gf2 operator- (const gf2 &other) const;
    gf2 &operator-= (const gf2 &other);

    gf2 operator* (const gf2 &other) const;
    gf2 &operator*= (const gf2 &other);
    gf2 operator/ (const gf2 &other) const;
    gf2 &operator/= (const gf2 &other);
    gf2 operator% (const gf2 &other) const;
    gf2 &operator%= (const gf2 &other);

    gf2 operator^ (int32_t degree) const;
    gf2 &operator^= (int32_t degree);
};

#endif