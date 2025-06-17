#include <cstdint>
#include <iostream>
#include "params.h"
#include "poly_mul.h"
#include "cpucycles.h"
#include "randombytes.h"
#include "sign.h"
#include "fips202.h"
#include "gf2.h"

void SSS_gf2_share(uint32_t s, const uint32_t *ids, uint32_t t, uint32_t n, uint32_t mod, uint32_t *shares)
{
    uint32_t poly[t];

    poly[0] = s;
    randombytes((unsigned char*)&poly[1], sizeof(uint32_t) * (t-1));

    for (uint32_t i = 0; i < n; i++) {
        shares[i] = 0;
        gf2 deg(mod, 1);
        for (uint32_t j = 0; j < t; j++) {
            shares[i] = (uint32_t)(gf2(mod, shares[i]) + gf2(mod, poly[j]) * deg);
            deg = deg * gf2(mod, ids[i]);
        }
    }
}

uint32_t SSS_gf2_recombine(uint32_t *shares, const uint32_t *ids, uint32_t t, uint32_t mod)
{
    uint32_t s = 0;
    
    for (uint32_t i = 0; i < t; i++) {
        gf2 l(mod, 1);
        for (uint32_t j = 0; j < t; j++) {
            if (j == i) {
                continue;
            }
            l *= gf2(mod, ids[j]) * ((gf2(mod, ids[j]) -  gf2(mod, ids[i])) ^ -1);
        }
        
        s = (uint32_t)(gf2(mod, s) + gf2(mod, shares[i]) * l);
    }
    
    return s;
}

void SSS_fp_share(uint32_t s, const uint32_t *ids, uint32_t t, uint32_t n, uint32_t p, uint32_t *shares)
{
    uint32_t poly[t];

    poly[0] = s;
    randombytes((unsigned char*)&poly[1], sizeof(uint32_t) * (t-1));

    for (uint32_t i = 0; i < t; i++) {
        poly[i] %= p;
        //std::cout << poly[i] << " ";
    }
    //std::cout << std::endl;

    for (uint32_t i = 0; i < n; i++) {
        shares[i] = 0;
        uint32_t deg = 1;
        for (uint32_t j = 0; j < t; j++) {
            shares[i] = ((uint64_t)shares[i] + (((uint64_t)poly[j] * deg) % p)) % p;
            deg = ((uint64_t)deg * ids[i]) % p;
        }
        //std::cout << shares[i] << " ";
    }
    //std::cout << std::endl;
}

int64_t gcdExtended(int64_t a, int64_t b, int64_t* x, int64_t* y)
{
    if (a == 0) {
        *x = 0, *y = 1;
        return b;
    }

    int64_t x1, y1;
    int64_t gcd = gcdExtended(b % a, a, &x1, &y1);

    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}

uint32_t inv_mod(int64_t a, uint32_t p)
{
    int64_t x, y;
    int64_t g = gcdExtended(a, p, &x, &y);
    if (abs(g) != 1){
        printf("Inverse doesn't exist g=%ld\n", g);
        return 0;
    }
    return ((g * x) % p + p) % p;
}

uint32_t SSS_fp_recombine(uint32_t *shares, const uint32_t *ids, uint32_t t, uint32_t p)
{
    uint32_t s = 0;
    
    for (uint32_t i = 0; i < t; i++) {
        uint32_t l = 1;
        for (uint32_t j = 0; j < t; j++) {
            if (j == i) {
                continue;
            }
            l = ((uint64_t)l * (((uint64_t)ids[j] * inv_mod((int64_t)ids[j] - ids[i], p)) % p)) % p;
        }
        //std::cout << l << " ";
        s = ((uint64_t)s + (((uint64_t)shares[i] * l) % p)) % p;
    }
    //std::cout << std::endl;
    return s;
}


void PQS_threshold_keygen_gf2(vk_t &vk, sk_t *sk, const uint32_t *ids, uint32_t tt, uint32_t n, uint32_t mod){
	uint32_t i,j;
	GenSeed(vk.seedA);
	
	polymatkl matA;
	GenMatrix(matA, vk.seedA);
	
    polyvecl s[n];
    for (uint32_t i = 0; i < n; i++) {
        GenSecret_s(s[i]);
    }
	
    polyveck t[n];
    for (uint32_t k = 0; k < n; k++) {
        matrixpoly(matA, s[k], t[k]);
        for(j=0; j<PQS_k; j++){
            for(i=0; i<PQS_n; i++){
                t[k].polynomial[j].coeffs[i] = round(t[k].polynomial[j].coeffs[i]);
            }
        }
    }
	
    for(j=0; j<PQS_k; j++){
        for(i=0; i<PQS_n; i++){
            vk.t.polynomial[j].coeffs[i] = 0;
            for (uint32_t k = 0; k < n; k++) {
                vk.t.polynomial[j].coeffs[i] += t[k].polynomial[j].coeffs[i];
                vk.t.polynomial[j].coeffs[i] &= (PQS_p-1);
            }
        }
    }

    sk_t sk_shares[n][n] = {0};
    for (uint32_t id1 = 0; id1 < n; id1++) {
        for(j=0; j<PQS_l; j++){
            for(i=0; i<PQS_n; i++){
                uint32_t shares[n] = {0}; 
                SSS_gf2_share(s[id1].polynomial[j].coeffs[i], ids, tt, n, mod, shares);
                for (uint32_t id2 = 0; id2 < n; id2++) {
                    sk_shares[id1][id2].s.polynomial[j].coeffs[i] = shares[id2];
                }
            }
        }
    }

    for (uint32_t id1 = 0; id1 < n; id1++) {
        for(j=0; j<PQS_l; j++){
            for(i=0; i<PQS_n; i++){
                sk[id1].s.polynomial[j].coeffs[i] = 0;
                for (uint32_t id2 = 0; id2 < n; id2++) {
                    sk[id1].s.polynomial[j].coeffs[i] = (uint32_t)(
                        gf2(mod, sk[id1].s.polynomial[j].coeffs[i]) 
                        + gf2(mod, sk_shares[id2][id1].s.polynomial[j].coeffs[i])
                    );
                }
            }
        }        
    }

    /*sk_t sk1;
	
    for(j=0; j<PQS_l; j++){
        for(i=0; i<PQS_n; i++){
            sk1.s.polynomial[j].coeffs[i] = 0;
            for (uint32_t k = 0; k < n; k++) {
                //sk1.s.polynomial[j].coeffs[i] += s[k].polynomial[j].coeffs[i];
                //sk1.s.polynomial[j].coeffs[i] &= (PQS_q-1);

                sk1.s.polynomial[j].coeffs[i] = (uint32_t)(
                    gf2(mod, sk1.s.polynomial[j].coeffs[i])
                    + gf2(mod, s[k].polynomial[j].coeffs[i])
                );
            }
        }
    }

    std::cout << "------------ CHECK :: VECTOR s1 ------------" << std::endl;
	for(i=0; i<PQS_l; i++){
		std::cout << "s[" << i << "] = ";
		print_polystruct(sk1.s.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;*/
	    
    /*std::cout << "------------ PUBLIC KEY :: VECTOR t -------------" << std::endl;
	for(i=0; i<PQS_k; i++){
		std::cout << "t[" << i << "] = ";
		print_polystruct(vk.t.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;

    polyveck t1;
    matrixpoly(matA, sk.s, t1);
	for(j=0; j<PQS_k; j++){
		for(i=0; i<PQS_n; i++){
			t1.polynomial[j].coeffs[i] = round(t1.polynomial[j].coeffs[i]);
		}
	}

    std::cout << "------------ CHECK :: VECTOR t1 -------------" << std::endl;
	for(i=0; i<PQS_k; i++){
		std::cout << "t1[" << i << "] = ";
		print_polystruct(t1.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;*/
}

void PQS_threshold_keygen_fp(vk_t &vk, sk_t *sk, const uint32_t *ids, uint32_t tt, uint32_t n, uint32_t p){
	uint32_t i,j;
	GenSeed(vk.seedA);
	
	polymatkl matA;
	GenMatrix(matA, vk.seedA);
	
    polyvecl s[n];
    for (uint32_t i = 0; i < n; i++) {
        GenSecret_s(s[i]);
    }
	
    polyveck t[n];
    for (uint32_t k = 0; k < n; k++) {
        matrixpoly(matA, s[k], t[k]);
        for(j=0; j<PQS_k; j++){
            for(i=0; i<PQS_n; i++){
                t[k].polynomial[j].coeffs[i] = round(t[k].polynomial[j].coeffs[i]);
            }
        }
    }
	
    for(j=0; j<PQS_k; j++){
        for(i=0; i<PQS_n; i++){
            vk.t.polynomial[j].coeffs[i] = 0;
            for (uint32_t k = 0; k < n; k++) {
                vk.t.polynomial[j].coeffs[i] += t[k].polynomial[j].coeffs[i];
                vk.t.polynomial[j].coeffs[i] &= (PQS_p-1);
            }
        }
    }

    sk_t sk_shares[n][n] = {0};
    for (uint32_t id1 = 0; id1 < n; id1++) {
        for(j=0; j<PQS_l; j++){
            for(i=0; i<PQS_n; i++){
                uint32_t shares[n] = {0}; 
                SSS_fp_share(s[id1].polynomial[j].coeffs[i], ids, tt, n, p, shares);
                for (uint32_t id2 = 0; id2 < n; id2++) {
                    sk_shares[id1][id2].s.polynomial[j].coeffs[i] = shares[id2];
                }
            }
        }
    }

    for (uint32_t id1 = 0; id1 < n; id1++) {
        for(j=0; j<PQS_l; j++){
            for(i=0; i<PQS_n; i++){
                sk[id1].s.polynomial[j].coeffs[i] = 0;
                for (uint32_t id2 = 0; id2 < n; id2++) {
                    sk[id1].s.polynomial[j].coeffs[i] = 
                        ((uint64_t)sk[id1].s.polynomial[j].coeffs[i] 
                        + sk_shares[id2][id1].s.polynomial[j].coeffs[i]) % p;
                }
                //sk[id1].s.polynomial[j].coeffs[i] &= (PQS_q-1);
            }
        }        
    }

    sk_t sk1;
	
    for(j=0; j<PQS_l; j++){
        for(i=0; i<PQS_n; i++){
            sk1.s.polynomial[j].coeffs[i] = 0;
            for (uint32_t k = 0; k < n; k++) {
                sk1.s.polynomial[j].coeffs[i] += s[k].polynomial[j].coeffs[i];
                sk1.s.polynomial[j].coeffs[i] &= (PQS_q-1);
            }
        }
    }

    std::cout << "------------ CHECK :: VECTOR s1 ------------" << std::endl;
	for(i=0; i<PQS_l; i++){
		std::cout << "s[" << i << "] = ";
		print_polystruct(sk1.s.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;
	    
    /*std::cout << "------------ PUBLIC KEY :: VECTOR t -------------" << std::endl;
	for(i=0; i<PQS_k; i++){
		std::cout << "t[" << i << "] = ";
		print_polystruct(vk.t.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;

    polyveck t1;
    matrixpoly(matA, sk.s, t1);
	for(j=0; j<PQS_k; j++){
		for(i=0; i<PQS_n; i++){
			t1.polynomial[j].coeffs[i] = round(t1.polynomial[j].coeffs[i]);
		}
	}

    std::cout << "------------ CHECK :: VECTOR t1 -------------" << std::endl;
	for(i=0; i<PQS_k; i++){
		std::cout << "t1[" << i << "] = ";
		print_polystruct(t1.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;*/
}

#define crh(OUT, IN, INBYTES) shake256(OUT, CRHBYTES, IN, INBYTES)

void PQS_threshold_sign(signat_t &sig, const unsigned char *m, uint32_t mlen, sk_t &sk, vk_t &vk){
	uint32_t i, j;
	unsigned char seed[SEEDBYTES];
	uint8_t mu[CRHBYTES];
	
	polymatkl matA;
	GenMatrix(matA, vk.seedA);	//sample matrix A
	
	uint32_t mod = (1 << (PQS_nu - PQS_mu));
	uint32_t delta1 = (1 << (PQS_nu-PQS_mu+1))*PQS_omega;
	uint32_t bound1 = (1 << (PQS_nu-PQS_d)) - 2*delta1 - 1;
	uint32_t bound2 = PQS_gamma - PQS_beta;
	uint32_t w1, z1;
	
	sgn:
	GenSeed(seed);
	// Generate intermediate vector y
	polyvecl y;
	GenSecret_y(y, seed);
	
	#ifdef DEBUG_MODE
	// Show vector y (l polynomials)
	cerr << "------------ SIGNING PROCEDURE :: VECTOR y -------------" << endl;
	for(int k=0; k<PQS_l; k++){
		cerr << "y[" << k << "] = ";
		print_polystruct(y.polynomial[k], PQS_n, PQS_q);
	}
	#endif
	
	// Prepare w=MSB(A*y, d)
	polyveck w;
	matrixpoly(matA, y, w);
	
	#ifdef DEBUG_MODE
	// Show vector w (k polynomials)
	cerr << "------------ SIGNING PROCEDURE :: VECTOR A*y -------------" << endl;
	for(int k=0; k<PQS_k; k++){
		cerr << "A*y[" << k << "] = ";
		print_polystruct(w.polynomial[k], PQS_n, PQS_q);
	}
	#endif
	
	for(j=0; j<PQS_k; j++)
		for(i=0; i<PQS_n; i++){
			w.polynomial[j].coeffs[i] = MSB(w.polynomial[j].coeffs[i], PQS_d);
		}
	
	#ifdef DEBUG_MODE
	// Show vector w (k polynomials)
	cerr << "------------ SIGNING PROCEDURE :: VECTOR MSB(A*y) -------------" << endl;
	for(int k=0; k<PQS_k; k++){
		cerr << "w[" << k << "] = ";
		print_polystruct(w.polynomial[k], PQS_n, PQS_q);
	}
	#endif
	
	// CRH(m) -> mu
	memset(mu, 0, sizeof(uint8_t)*CRHBYTES);
	crh(mu, m, mlen);
	#ifdef DEBUG_MODE
	// Show a hash mu
	cerr << "------------ SIGNING PROCEDURE :: CRH(Message) -------------" << endl;
	for(int k=0; k<CRHBYTES; k++){
		cerr << mu[k] << " ";
	}
	cerr << endl;
	#endif	
	
	// Generate polynomial c
	challenge(sig.c, &mu[0], w);
	
	#ifdef DEBUG_MODE
	// Show polynomial c (challenge)
	cerr << "------------ SIGNING PROCEDURE :: POLYNOMIAL c -------------" << endl;
		print_polystruct(sig.c, PQS_n, PQS_q);
	#endif
	
	// Calculate z
	for(i=0; i<PQS_l; i++){
		POLYMUL_MODE(sk.s.polynomial[i], sig.c, sig.z.polynomial[i], PQS_n, PQS_q);
		// add y[i] to z[i]
		for(j=0; j<PQS_n; j++){
			sig.z.polynomial[i].coeffs[j] += y.polynomial[i].coeffs[j];
			sig.z.polynomial[i].coeffs[j] &= (PQS_q - 1); // Reduce modulo q
		}
	}
	
	#ifdef DEBUG_MODE
	// Show vector z (l polynomials)
	cerr << "------------ SIGNING PROCEDURE :: VECTOR z -------------" << endl;
	for(int k=0; k<PQS_l; k++){
		cerr << "z[" << k << "] = ";
		print_polystruct(sig.z.polynomial[k], PQS_n, PQS_q);
	}
	#endif
	
	// Calculate w (memory is already allocated, we don't need old w anymore
	matrixpoly(matA, sig.z, w);
	
	poly subpoly;
	for(i=0; i<PQS_k; i++){
		POLYMUL_MODE(vk.t.polynomial[i], sig.c, subpoly, PQS_n, PQS_q);
		for(j=0; j<PQS_n; j++){
			subpoly.coeffs[j] *= mod;
			w.polynomial[i].coeffs[j] -= subpoly.coeffs[j];
			w.polynomial[i].coeffs[j] &= (PQS_q-1);
		}
	}
	
	#ifdef DEBUG_MODE
	// Show vector w (k polynomials)
	cerr << "------------ SIGNING PROCEDURE :: VECTOR w -------------" << endl;
	for(int k=0; k<PQS_k; k++){
		cerr << "w[" << k << "] = ";
		print_polystruct(w.polynomial[k], PQS_n, PQS_q);
	}
	#endif
	
	for(i=0; i<PQS_k; i++){ // Check w norm bound
		for(j=0; j<PQS_n; j++){
			w1 = LSB(w.polynomial[i].coeffs[j], PQS_nu-PQS_d);
			w1 -= (delta1 + 1);
			if(w1 >= bound1){
				//cerr << "PQS_sign() RESTARTING (bad w) ... " << endl;
				goto sgn;
			}
		}
	}
	
	for(i=0; i<PQS_l; i++){ // Check z norm bound
		for(j=0; j<PQS_n; j++){
			z1 = sig.z.polynomial[i].coeffs[j] - bound2;
			if(z1 <= PQS_q - 2*bound2 - 1){
				//cerr << "PQS_sign() RESTARTING (bad z) ... " << endl;
				goto sgn;
			}
		}
	}
}

void test_SSS()
{
    uint32_t s = 123; 
    uint32_t p = 4294967291;

    uint32_t ids[5] = {1,2,3,4,5};
    uint32_t shares[5] = {0};

    SSS_fp_share(s, ids, 3, 5, p, shares);
    
    uint32_t ids1[3] = {ids[0], ids[4], ids[1]};
    uint32_t shares1[3] = {shares[0], shares[4], shares[1]};

    std::cout << SSS_fp_recombine(shares1, ids1, 3, p) << std::endl;
}

void test_PQS_threshold_keygen()
{
    uint32_t t = 3;
    uint32_t n = 5;
    uint32_t p = 4294967291;

    vk_t vk;
    sk_t sk_shares[n];
    uint32_t ids[n] = {0};

    for (uint32_t i = 0; i < n; i++) {
        ids[i] = i + 1;
    }

    PQS_threshold_keygen_fp(vk, sk_shares, ids, t, n, p);
    
    uint32_t shares_ind[t] = {0, 4, 1};
    uint32_t ids1[t];
    for (uint32_t i = 0; i < t; i++) {
        ids1[i] = ids[shares_ind[i]];
    }

    sk_t sk_shares1[t];
    for (uint32_t i = 0; i < t; i++) {
        sk_shares1[i] = sk_shares[shares_ind[i]];
    }

    sk_t sk;

    for(uint32_t j=0; j<PQS_l; j++){
        for(uint32_t i=0; i<PQS_n; i++){
            uint32_t shares[t];
            for (uint32_t k = 0; k < t; k++) {
                shares[k] = sk_shares1[k].s.polynomial[j].coeffs[i];
            }
            sk.s.polynomial[j].coeffs[i] = SSS_fp_recombine(shares, ids1, t, p) & (PQS_q-1);
        }
    }
    
    std::cout << "------------ SECRET KEY :: VECTOR s ------------" << std::endl;
	for(uint32_t i=0; i<PQS_l; i++){
		std::cout << "s[" << i << "] = ";
		print_polystruct(sk.s.polynomial[i], PQS_n, PQS_q);
	}
    std::cout << std::endl;

    unsigned char m[] = "My test message";
    signat_t sig;
    PQS_sign(sig, &m[0], sizeof(m), sk, vk);
    
    bool success = PQS_verify(sig, &m[0], sizeof(m), vk);
    if(success){
        std::cout << "Verify Ok" << std::endl;
    } else {
        std::cout << "Verify fail" << std::endl;
    }
}

void test_gf2() 
{
/*
Код для sage для поиска примитивых полиномов

i = (1<<23)+1
while(i < 1<<24):
    pol = 1
    for j in range(1,24):
        pol += ((i>>j)&1) * x^j
    if pol.is_primitive():
        print(i, pol)
    i += 2

*/

    gf2 a(7, 2), b(7, 3);

    std::cout << (uint32_t)a << " + " << (uint32_t)b << " = " << (uint32_t)(a+b) << std::endl;

    uint32_t mod = (1<<5) + (1<<2) + 1;
    gf2 g(mod, 3);
    std::cout << (uint32_t)g << ", ";
    for (gf2 x = g*g; x != g && x != gf2(mod, 0); x *= g){
        std::cout << (uint32_t)x << ", ";
    }
    std::cout << std::endl;

    gf2 c(mod, 20), d(mod, 3);
    std::cout << (uint32_t)c << " / " << (uint32_t)d << " = " << (uint32_t)(c/d) << std::endl;

    gf2 e(mod, 3);
    int32_t deg = -2;
    std::cout << (uint32_t)e << "^" << deg << " = " << (uint32_t)(e^deg) << std::endl;
    
    uint32_t mod23 = (1<<23) + (1<<14) + 1;
    gf2 v(mod23, 1234567);
    std::cout << (uint32_t)v << "*" << (uint32_t)(v^-1) << " = " << (uint32_t)(v * (v^-1)) << std::endl;
}

int main()
{
    test_PQS_threshold_keygen();

    return 0;
}