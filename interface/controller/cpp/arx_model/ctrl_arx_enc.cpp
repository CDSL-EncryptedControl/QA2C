// get tcp_protocol decription
#include "tcp_protocol_client.h" 
// #include "tcp_protocol_client_windows.h" 

// get other tools
#include <iostream>
// #include <cstdint> //if use Windows environment, then active this.
#include <string>
#include <chrono>
#include <cmath>
using namespace std;

// get model description
#include "model_enc.h"

// init tcp host and port
const string host = "172.25.32.1";
const int port = 9999;

int main()
{
    // set simulation(this section have to set same with plant)
    double sampling_time = 0.05;
    bool run_signal = true;

    // get crypto model from model_enc.h
    crypto crypto_cl = crypto();
    enc_for_arx enc_4_arx = enc_for_arx(crypto_cl);
    enc_4_arx.set_level(2048, 2048);
    arx_enc arx_enc_v = arx_enc(crypto_cl.get_crypto(), crypto_cl.get_relinkey(), crypto_cl.get_galoiskeys());
    arx_enc_v.set_hyper_params(enc_4_arx.get_HG_enc(), enc_4_arx.get_HL_enc(), enc_4_arx.get_Href_enc());
    arx_enc_v.set_io(enc_4_arx.get_y_enc(), enc_4_arx.get_u_enc(), enc_4_arx.get_ref_enc());

    // set tcp client
    tcp_client tccp = tcp_client(host, port);
    string signal;

    // for check cycle time
    auto stc = chrono::high_resolution_clock::now();
    auto edc = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    double run_time;

    // input/output initialization
    vector<double> y(2, 0.0);
    vector<double> u(2, 0.0);
    vector<int64_t> int_u(2, 0LL);

    // refernece setting
    vector<double> ref(2, 0.0);
    vector<double> ref_delta(2, 0.0);
    ref[0] = 15 * M_PI / 180;
    ref[1] = -15 * M_PI / 180;

    while(run_signal)
    {
        signal = tccp.Recv<string>();

        if(signal == "run")
        {
            // start clock set
            stc = chrono::high_resolution_clock::now();

            // get plant output
            double y0 = tccp.Recv<double>();
            double y1 = tccp.Recv<double>();
            tccp.Recv<double>();
            tccp.Recv<double>();
            y[0] = y0;
            y[1] = y1;
            
            // send control input data
            tccp.Send<double>(u[0]);
            tccp.Send<double>(u[1]);

            // send ref
            tccp.Send<double>(ref[0]);
            tccp.Send<double>(ref[1]);
            ref_delta[0] = ref[0] * sampling_time;
            ref_delta[1] = ref[1] * sampling_time;

            // y and u value encryption after packing
            Ciphertext signal_y = enc_4_arx.enc_signal_4_y(y);
            Ciphertext signal_u = enc_4_arx.enc_signal_4_u(u);
            Ciphertext signal_ref = enc_4_arx.enc_signal_4_ref(ref_delta);

            // -- controller description -- //
            // ========================================================== //
            // ctrl mem update on encrypted space after encryption input/output value
            arx_enc_v.mem_update(signal_y, signal_u, signal_ref);

            // get control input on ciphertext space
            Ciphertext enc_u = arx_enc_v.get_output();
            // ========================================================== //

            int_u = enc_4_arx.dec_signal(enc_u); 

            u[0] = (double)(int_u[0]) / enc_4_arx.get_level()[0] / enc_4_arx.get_level()[1];
            u[1] = (double)(int_u[4]) / enc_4_arx.get_level()[0] / enc_4_arx.get_level()[1];

            // end clock set
            edc = chrono::high_resolution_clock::now();
            duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
            run_time = duration.count() / 1000000;
            cout << "loop time: " << run_time << "ms" << endl;

        }
        else if(signal == "end")
        {
            // end of loop signal get
            run_signal = false;
            break;
        }
    }

    return 0;
}
