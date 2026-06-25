#include "seal/seal.h"
// #include <cstdint> //if use Windows environment, then active this.
#include <vector>
using namespace seal;
using namespace std;

class crypto
{
    private:
        // cryptocontext for encryption
        EncryptionParameters parms;
        SEALContext context;

        // key_pair for encryption
        KeyGenerator keygen; 
        SecretKey secret_key;
        PublicKey public_key;
        RelinKeys relin_keys;
        GaloisKeys galois_keys;
        Encryptor encryptor;
        Decryptor decryptor;
        Evaluator evaluator;
        BatchEncoder batch_encoder;
        size_t slot_count;

        // parameters helper (parameters setting)
        EncryptionParameters create_paramters()
        {
            EncryptionParameters parms(scheme_type::bgv);
            size_t poly_modulus_degree = 8192;
            parms.set_poly_modulus_degree(poly_modulus_degree);
            parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {57, 57, 57}));
            parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 42));
            
            return parms;
        };

        PublicKey create_public_key(KeyGenerator& keygen)
        {
            PublicKey publickey;
            keygen.create_public_key(publickey);
            
            return publickey;
        };

    public:
        crypto(): 
            parms(create_paramters()),
            context(this->parms),
            keygen(this->context),
            secret_key(this->keygen.secret_key()),
            public_key(create_public_key(this->keygen)),
            decryptor(this->context, this->secret_key),
            encryptor(this->context, this->public_key),
            evaluator(this->context),
            batch_encoder(this->context)
        {
            this->keygen.create_relin_keys(this->relin_keys);
            
            vector<int> steps = {1, 2};
            this->keygen.create_galois_keys(steps, this->galois_keys);

            this->slot_count = this->batch_encoder.slot_count();
        };

        ~crypto(){};

        SEALContext get_crypto()
        {
            return this->context;
        };

        RelinKeys get_relinkey()
        {
            return this->relin_keys;
        };

        GaloisKeys get_galoiskeys()
        {
            return this->galois_keys;
        };


        size_t get_slotsize()
        {
            return this->slot_count;
        };

        Ciphertext enc_vector(vector<int64_t> vec)
        {
            // make plaintext with packing
            Plaintext plain;
            this->batch_encoder.encode(vec, plain);

            // make ciphertext with encryptor
            Ciphertext cipher;
            this->encryptor.encrypt(plain, cipher);

            return cipher;
        };

        vector<int64_t> dec_ciphertext(Ciphertext cipher)
        {
            // make plaintext with decryptor
            Plaintext plain;
            this->decryptor.decrypt(cipher, plain);

            // make vector with unpacking
            vector<int64_t> vec;
            this->batch_encoder.decode(plain, vec);

            return vec;
        };
};

class enc_for_arx
{
    private:
        crypto& crypto_cl;

        int64_t r = 1000;
        int64_t s = 1000;

        // you can get HG_q and HL_q from controller/py/controller.py file
        int64_t HG_q[12][2] = {{-47,-9},
                            {-19,-3},
                            {-1247,-239},
                            {-528,-101},
                            {-32837,-5857},
                            {-14062,-2964},
                            {-669932,495004},
                            {-482074,-439352},
                            {1513507,-1026680},
                            {1063061,951314},
                            {-810642,538650},
                            {-567310,-509916}};
        int64_t HL_q[12][2] = {{0,0},
                            {0,0},
                            {0,0},
                            {0,0},
                            {2,7},
                            {0,3},
                            {179,-62},
                            {-39,241},
                            {-1007,1},
                            {-72,-855},
                            {2874,53},
                            {111,2658}};
        int64_t Href[12][6] = {{-49,-9,2,0,0,0},
                            {-20,-4,1,0,0,0},
                            {-1304,-252,64,12,3,1},
                            {-553,-107,27,5,1,0},
                            {-34363,-6185,1693,311,89,27},
                            {-14711,-3112,724,154,38,14},
                            {-706766,498673,35624,-21856,1141,-3250},
                            {-501202,-450119,24252,20806,1697,2595},
                            {931523,-644716,44486,-27578,-40403,35113},
                            {666161,622731,28973,20754,-33406,-38463},
                            {-189037,152490,-81870,49109,63182,-49249},
                            {-149672,-169387,-53978,-41721,50352,56363}};

        // encrypted gain that packed and encrypted like PQ_enc[0] = {HG_q[0, 0], HG_q[0, 1], HL_q[0, 0]}                              
        vector<Ciphertext> HG_enc;
        vector<Ciphertext> HL_enc;
        vector<Ciphertext> Href_enc;

        // encrypted signal sequence that packed and encrypted like S_enc[0] = {y[0, 0], y[1, 0], u[0, 0]}
        Ciphertext y_signal_enc;
        Ciphertext u_signal_enc;
        Ciphertext ref_signal_enc;
        vector<Ciphertext> y_enc;
        vector<Ciphertext> u_enc;
        vector<Ciphertext> ref_enc;
    
    public:
        enc_for_arx(crypto& crypto_class): 
            crypto_cl(crypto_class),
            HG_enc(6),
            HL_enc(6),
            Href_enc(6),
            y_enc(6),
            u_enc(6),
            ref_enc(1)
        {
            for(int i = 0; i < 6; i++)
            {
                vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);

                vec[0] = this->HG_q[2*i][0];
                vec[1] = this->HG_q[2*i][1];
                vec[4] = this->HG_q[2*i+1][0];
                vec[5] = this->HG_q[2*i+1][1];

                this->HG_enc[i] = this->crypto_cl.enc_vector(vec);
            }

            for(int i = 0; i < 6; i++)
            {
                vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);

                vec[0] = this->HL_q[2*i][0];
                vec[1] = this->HL_q[2*i][1];
                vec[4] = this->HL_q[2*i+1][0];
                vec[5] = this->HL_q[2*i+1][1];

                this->HL_enc[i] = this->crypto_cl.enc_vector(vec);
            }

            for(int i = 0; i < 6; i++)
            {
                vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);

                vec[0] = this->Href[2*i][4];
                vec[1] = this->Href[2*i][5];
                vec[4] = this->Href[2*i+1][4];
                vec[5] = this->Href[2*i+1][5];

                this->Href_enc[i] = this->crypto_cl.enc_vector(vec);
            }

            for(int i = 0; i < 6; i++)
            {
                vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);

                this->y_enc[i] = this->crypto_cl.enc_vector(vec);
            }

            for(int i = 0; i < 6; i++)
            {
                vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);

                this->u_enc[i] = this->crypto_cl.enc_vector(vec);
            }

            vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);
            this->ref_enc[0] = this->crypto_cl.enc_vector(vec);
        };

        ~enc_for_arx(){};

        void set_level(int64_t r, int64_t s)
        {
            this->r = r;
            this->s = s;
        };

        vector<int64_t> get_level()
        {
            vector<int64_t> level(2);
            level[0] = this->r;
            level[1] = this->s;
            
            return level;
        };

        vector<Ciphertext> get_HG_enc()
        {
            return this->HG_enc;
        };

        vector<Ciphertext> get_HL_enc()
        {
            return this->HL_enc;
        };

        vector<Ciphertext> get_Href_enc()
        {
            return this->Href_enc;
        };

        vector<Ciphertext> get_y_enc()
        {
            return this->y_enc;
        };

        vector<Ciphertext> get_u_enc()
        {
            return this->u_enc;
        };

        vector<Ciphertext> get_ref_enc()
        {
            return this->ref_enc;
        };

        Ciphertext enc_signal_4_y(vector<double> y_q)
        {
            vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);
            vec[0] = (int)(y_q[0] * this->r);
            vec[1] = (int)(y_q[1] * this->r);
            vec[4] = (int)(y_q[0] * this->r);
            vec[5] = (int)(y_q[1] * this->r);

            this->y_signal_enc = this->crypto_cl.enc_vector(vec);

            return this->y_signal_enc;
        };

        Ciphertext enc_signal_4_u(vector<double> u_q)
        {
            vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);
            vec[0] = (int)(u_q[0] * this->r);
            vec[1] = (int)(u_q[1] * this->r);
            vec[4] = (int)(u_q[0] * this->r);
            vec[5] = (int)(u_q[1] * this->r);

            this->u_signal_enc = this->crypto_cl.enc_vector(vec);

            return this->u_signal_enc;
        };

        Ciphertext enc_signal_4_ref(vector<double> ref_q)
        {
            vector<int64_t> vec(this->crypto_cl.get_slotsize(), 0LL);
            vec[0] = (int)(ref_q[0] * this->r);
            vec[1] = (int)(ref_q[1] * this->r);
            vec[4] = (int)(ref_q[0] * this->r);
            vec[5] = (int)(ref_q[1] * this->r);

            this->y_signal_enc = this->crypto_cl.enc_vector(vec);

            return this->y_signal_enc;
        };

        vector<int64_t> dec_signal(Ciphertext cipher)
        {
            vector<int64_t> vec;
            vec = this->crypto_cl.dec_ciphertext(cipher);

            return vec;
        }
};

class arx_enc
{
    private:
        SEALContext ctext;

        Evaluator evaluator;

        RelinKeys relin_keys;

        GaloisKeys galois_keys;

        vector<Ciphertext>hg;
        vector<Ciphertext>hl;
        vector<Ciphertext>href;

        vector<Ciphertext>io_y; // oldest -> newest
        vector<Ciphertext>io_u;
        vector<Ciphertext>io_ref;

    public:
        arx_enc(const SEALContext& context, RelinKeys relin_keys, GaloisKeys galois_keys): 
            ctext(context), 
            evaluator(context), 
            hg(6), 
            hl(6),
            href(6),
            io_y(6), 
            io_u(6),
            io_ref(1)
        {
            this->relin_keys = relin_keys;
            this->galois_keys = galois_keys;
        };
        ~arx_enc(){};

        void set_hyper_params(const vector<Ciphertext>& hg, const vector<Ciphertext>& hl, const vector<Ciphertext>& href)
        {
            for(int i = 0; i < 6; i++)
            {
                this->hg[i] = hg[i];
                this->hl[i] = hl[i];
                this->href[i] = href[i];
            }
        }

        void set_io(const vector<Ciphertext>& io_y, const vector<Ciphertext>& io_u, const vector<Ciphertext>& io_ref)
        {
            for(int i = 0; i < 6; i++)
            {
                this->io_y[i] = io_y[i];
                this->io_u[i] = io_u[i];
            }
            this->io_ref[0] = io_ref[0];
        }

        void mem_update(Ciphertext new_y, Ciphertext new_u, Ciphertext new_ref)
        {
            io_y.erase(io_y.begin());
            io_y.push_back(new_y);

            io_u.erase(io_u.begin());
            io_u.push_back(new_u);

            io_ref.erase(io_ref.begin());
            io_ref.push_back(new_ref);
        };

        Ciphertext get_output()
        {
            vector<Ciphertext> hg_mul(6);
            vector<Ciphertext> hl_mul(6);
            vector<Ciphertext> href_mul(6);
            Ciphertext t_mul;
            Ciphertext r_sum;
            Ciphertext u_enc;

            for(int i = 0; i < 6; i++)
            {
                this->evaluator.multiply(this->hg[i], this->io_y[i], t_mul);
                this->evaluator.relinearize_inplace(t_mul, this->relin_keys);
                hg_mul[i] = t_mul;

                this->evaluator.multiply(this->hl[i], this->io_u[i], t_mul);
                this->evaluator.relinearize_inplace(t_mul, this->relin_keys);
                hl_mul[i] = t_mul;

                this->evaluator.multiply(this->href[i], this->io_ref[0], t_mul);
                this->evaluator.relinearize_inplace(t_mul, this->relin_keys);
                href_mul[i] = t_mul;
            }

            r_sum = hg_mul[0];

            for(int i = 1; i < 18; i++)
            {
                if(i < 6)
                {
                    this->evaluator.add_inplace(r_sum, hg_mul[i]);
                }
                else if(i < 12)
                {
                    this->evaluator.add_inplace(r_sum, hl_mul[i-6]);
                }
                else
                {
                    this->evaluator.add_inplace(r_sum, href_mul[i-12]);
                }
            }

            u_enc = r_sum;
            this->evaluator.rotate_rows_inplace(r_sum, 1, this->galois_keys);
            this->evaluator.add_inplace(u_enc, r_sum);

            return u_enc;
        };
};
