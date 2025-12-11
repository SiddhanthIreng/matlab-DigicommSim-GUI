# Digital Communication System Simulator

This project is a MATLAB-based graphical user interface (GUI) application designed to simulate and visualize the core stages of a digital communication link, including Source Coding, Channel Coding, Modulation, Channel Effects (AWGN), and Decoding.The simulator is an excellent tool for students and engineers to observe the trade-offs between different communication parameters (e.g., SNR, modulation schemes, and error correction) in real-time.
# Features

**Audio Processing:** Quantizes and processes standard .wav audio files or audio recording input.<br/>
**Quantization:** Supports 4-bit and 8-bit Pulse Code Modulation (PCM).<br/>
**Source Coding:** Includes an option for Huffman (Variable Length) coding.<br/>
**Channel Coding:** Implements Hamming (7,4) error correction coding.<br/>
**Modulation Schemes:** Supports BPSK, QPSK, 16-QAM, and 64-QAM.<br/>
**Channel Effects:** Simulates the effects of Additive White Gaussian Noise (AWGN) with adjustable Signal-to-Noise Ratio (SNR).<br/>
**Visualization:** Real-time plots of Transmitted/Received Waves, Constellation Diagrams, and Coding/Decoding Tables.<br/>
# Getting Started

**Prerequisites**

This project requires MATLAB to run.<br/>

**Installation (Download)**

Clone the repository:  

    git clone https://github.com/SiddhanthIreng/matlab-DigicommSim-GUI.git
Navigate to the project directory:  

    cd matlab-DigicommSim-GUI
Download as ZIP: Alternatively, you can download the repository as a ZIP file and extract the contents to a folder on your computer.

# How to Use the Simulator
The entire application runs from a single main script/app file (likely named DigitalCommSystem.m or similar).

**Launch the Application**

1. Open MATLAB.
2. Navigate to the matlab-DigicommSim-GUI folder in the MATLAB Current Folder pane.
3. Type the name of the main file (e.g., DigitalCommSystem) and press Enter, or simply double-click the main app file (if it's a .m file).
   
 Or, 
 if you are using MATLAB from browser.
 
  1. Log in to https://matlab.mathworks.com
  2. Upload the file:
    Go to Matlab drive and click Upload button.<br/>Select DigitalCommSystem.m file and the audio .wav files(optional) from the download folder and run.

*Note: Dont change the .m file name as the class name is same.*

**Simulation**

  1. Load Source Audio: Click the "Load Audio" button and select a .wav file or cick the "Record Audio" button to record, set the recording duration(upto 5 sec).
  2. Configure Parameters: <br/>Adjust the Quantization (4-bit or 8-bit).<br/>Select the Modulation Scheme (BPSK, QPSK, etc.).<br/>Enable/Disable Source (Huffman) and Channel (Hamming) coding.
  3. Set Channel Conditions: Use the SNR Slider to set the desired Signal-to-Noise Ratio for the AWGN channel.
  4. Run Simulation: Click the "RUN SIMULATION" button.
  5. Analyze Results: View the different tabs to analyze result
  6. Listen: Use the "Play Output" button to hear the decoded audio and evaluate the quality loss due to the channel and quantization, and use the "Play Input" to hear your original Input audio.
# License
  This project is licensed under the MIT License - see the LICENSE file for details.
