classdef DigitalCommSystem < handle
    % DIGITALCOMMSYSTEM - End-to-End Digital Communication Simulator 
    properties
        % GUI Components
        Fig
        Grid
        TabGroup
        
        % Panels
        ControlPanel
        
        % Widgets
        RecordBtn
        PlayOriginalBtn
        PlayRxBtn
        LoadBtn
        DurationField
        
        % Settings
        QuantizationDropDown
        ModulationDropDown
        chkSource   % Toggle for Huffman
        chkChannel  % Toggle for Hamming
        
        SnrSlider
        SnrLabel
        StatusLabel
        RunBtn
        
        % Data
        OriginalAudio
        Fs = 8000;
        RecObj
        ProcessedAudio
        
        % Simulation Data
        HuffDict        % Huffman Dictionary
        TxSymbols       % Raw integer symbols (before source coding)
        TxBits_Src      % Bits after source coding
        TxBits_Ch       % Bits after channel coding
        RxBits_Ch       % Received bits
        RxBits_Src      % Bits after channel decoding
        RxSymbols       % Recovered integer symbols
        TxSym
        ModType
        QLevel
        HammingPadding  % Number of zeros added for (7,4) alignment
        
        % Tabs
        TabSource
        TabADC
        TabCoding
        TabModulation
        TabDemod
        TabDecoding
        TabOutput
    end
    
    methods
        function app = DigitalCommSystem()
            app.createGUI();
        end
        
        function createGUI(app)
            app.Fig = uifigure('Name', 'Digital Comm System', 'Position', [50,50,1300,800]);
            app.Grid = uigridlayout(app.Fig, [1,2]);
            app.Grid.ColumnWidth = {280, '1x'};
            
            app.ControlPanel = uipanel(app.Grid, 'Title', 'Control Center');
            app.ControlPanel.Layout.Row = 1; app.ControlPanel.Layout.Column = 1;
            ctrlGrid = uigridlayout(app.ControlPanel, [16,1]);
            ctrlGrid.RowHeight = {30,30,30,20,30,30,30,30,30,30,30,20,30,30,30,'1x'};
            
            % Source settings
            uilabel(ctrlGrid,'Text','1. Source Settings','FontWeight','bold','BackgroundColor',[0.9 0.9 0.9]);
            g1 = uigridlayout(ctrlGrid,[1,2]); g1.ColumnWidth = {'1x',50}; g1.Padding=[0 0 0 0];
            uilabel(g1,'Text','Duration (s):'); app.DurationField = uieditfield(g1,'numeric','Value',2,'Limits',[1 5]);
            app.RecordBtn = uibutton(ctrlGrid,'Text','Record Audio','ButtonPushedFcn',@app.onRecord);
            app.LoadBtn = uibutton(ctrlGrid,'Text','Load .WAV File','ButtonPushedFcn',@app.onLoad);
            
            uilabel(ctrlGrid,'Text','-------------------------');
            uilabel(ctrlGrid,'Text','2. Coding & Modulation','FontWeight','bold','BackgroundColor',[0.9 0.9 0.9]);
            app.chkSource = uicheckbox(ctrlGrid,'Text','Enable Huffman Coding','Value',true);
            app.chkChannel = uicheckbox(ctrlGrid,'Text','Enable Hamming (7,4)','Value',true);
            
            uilabel(ctrlGrid,'Text','Quantization Level:');
            app.QuantizationDropDown = uidropdown(ctrlGrid,'Items',{'4-bit (16 levels)','8-bit (256 levels)'},'Value','8-bit (256 levels)');
            
            uilabel(ctrlGrid,'Text','Modulation Scheme:');
            app.ModulationDropDown = uidropdown(ctrlGrid,'Items',{'BPSK','QPSK','16-QAM','BFSK','ASK'},'Value','BPSK');
            
            g2 = uigridlayout(ctrlGrid, [1,2]); g2.ColumnWidth = {'1x',60}; g2.Padding = [0 0 0 0];
            app.SnrLabel = uilabel(g2, 'Text', 'SNR: 15 dB');
            app.SnrSlider = uislider(ctrlGrid, 'Limits', [-5,30], 'Value', 15, 'ValueChangedFcn', @(s,~)app.updateSnrLabel(s));
            
            uilabel(ctrlGrid,'Text','-------------------------');
            app.RunBtn = uibutton(ctrlGrid,'Text','RUN SIMULATION','FontWeight','bold','BackgroundColor',[0.3 0.6 1],'FontColor','white','ButtonPushedFcn',@app.runSimulation);
            app.StatusLabel = uilabel(ctrlGrid,'Text','Status: Ready');
            
            g3 = uigridlayout(ctrlGrid,[1,2]);
            app.PlayOriginalBtn = uibutton(g3,'Text','Play Input','Enable','off','ButtonPushedFcn',@app.playOriginal);
            app.PlayRxBtn = uibutton(g3,'Text','Play Output','Enable','off','ButtonPushedFcn',@app.playProcessed);
            
            % Right tabs
            app.TabGroup = uitabgroup(app.Grid); app.TabGroup.Layout.Row = 1; app.TabGroup.Layout.Column = 2;
            app.TabSource = uitab(app.TabGroup,'Title','1. Source');
            app.TabADC = uitab(app.TabGroup,'Title','2. ADC & Bits');
            app.TabCoding = uitab(app.TabGroup,'Title','3. Coding Tables');
            app.TabModulation = uitab(app.TabGroup,'Title','4. Mod Waves');
            app.TabDemod = uitab(app.TabGroup,'Title','5. Demodulation');
            app.TabDecoding = uitab(app.TabGroup,'Title','6. Decoding Tables');
            app.TabOutput = uitab(app.TabGroup,'Title','7. Comparison');
            
            % Default tone
            t = 0:1/app.Fs:0.5; app.OriginalAudio = sin(2*pi*440*t)';
        end
        
        function updateSnrLabel(app, src)
            app.SnrLabel.Text = sprintf('SNR: %.1f dB', src.Value);
        end
        
        function onRecord(app, ~, ~)
            duration = app.DurationField.Value;
            app.StatusLabel.Text = 'Recording...'; drawnow;
            try
                app.RecObj = audiorecorder(app.Fs, 16, 1);
                recordblocking(app.RecObj, duration);
                app.OriginalAudio = getaudiodata(app.RecObj);
                app.StatusLabel.Text = 'Recorded.'; app.PlayOriginalBtn.Enable = 'on'; app.plotSource();
            catch
                app.StatusLabel.Text = 'Mic Error';
            end
        end
        
        function onLoad(app, ~, ~)
            [file,path] = uigetfile('*.wav');
            if isequal(file,0), return; end
            [y,fs] = audioread(fullfile(path,file));
            if fs~=app.Fs, y = resample(y, app.Fs, fs); end
            if size(y,2)>1, y = mean(y,2); end
            if length(y) > 5*app.Fs, y = y(1:5*app.Fs); end
            app.OriginalAudio = y; app.PlayOriginalBtn.Enable='on'; app.StatusLabel.Text='Loaded.'; app.plotSource();
        end
        
        function runSimulation(app, ~, ~)
            if isempty(app.OriginalAudio)
                uialert(app.Fig, 'Load Audio first', 'Error'); return;
            end
            app.StatusLabel.Text = 'Processing...'; drawnow;
            try
                % ---- Quantization ----
                max_amp = max(abs(app.OriginalAudio));
                if max_amp > 1e-6, x = app.OriginalAudio / max_amp; else x = app.OriginalAudio; end
                if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), qBit = 4; else qBit = 8; end
                L = 2^qBit;
                x_int = floor(((x + 1)/2)*(L-1));
                x_int(x_int<0)=0; x_int(x_int> L-1)=L-1;
                symbols = double(x_int);
                raw_bin_strs = dec2bin(x_int, qBit);
                
                app.TxSymbols = symbols; % Store for visualization
                
                % ---- Source coding (Huffman) ----
                if app.chkSource.Value
                    tabS = tabulate(symbols);
                    unique_sym = tabS(tabS(:,2)>0,1);
                    probs = tabS(tabS(:,2)>0,3)/100;
                    app.HuffDict = huffmandict(unique_sym, probs);
                    app.TxBits_Src = huffmanenco(symbols, app.HuffDict);
                else
                    app.HuffDict = {};
                    app.TxBits_Src = reshape(de2bi(symbols, qBit, 'left-msb').', [], 1);
                end
                app.TxBits_Src = double(app.TxBits_Src(:));
                
                % ---- Channel coding (Hamming) ----
                if app.chkChannel.Value
                    n=7; k=4;
                    bits_in = app.TxBits_Src;
                    remBits = mod(length(bits_in), k); app.HammingPadding = 0;
                    if remBits>0, app.HammingPadding = k-remBits; bits_in = [bits_in; zeros(app.HammingPadding,1)]; end
                    msg_words = reshape(bits_in, k, []).';
                    coded_words = encode(msg_words, n, k, 'hamming/binary');
                    app.TxBits_Ch = reshape(coded_words.', [], 1);
                else
                    app.HammingPadding = 0;
                    app.TxBits_Ch = app.TxBits_Src;
                end
                app.TxBits_Ch = double(app.TxBits_Ch(:));
                
                % ---- Modulation ----
                app.ModType = app.ModulationDropDown.Value;
                M = 2; 
                switch app.ModType
                    case 'BPSK', M=2;
                    case 'QPSK', M=4;
                    case '16-QAM', M=16;
                    case 'BFSK', M=2;
                    case 'ASK', M=2; % OOK (On-Off Keying)
                end
                bitsIn = app.TxBits_Ch;
                bps = log2(M);
                modPad = 0;
                if mod(length(bitsIn), bps) ~= 0, modPad = bps - mod(length(bitsIn), bps); bitsIn = [bitsIn; zeros(modPad,1)]; end
                
                switch app.ModType
                    case 'BPSK'
                        tx_sym = pskmod(bitsIn, M);
                    case 'QPSK'
                        tx_sym = pskmod(bitsIn, M, pi/4, 'InputType', 'bit');
                    case '16-QAM'
                        tx_sym = qammod(bitsIn, M, 'InputType', 'bit', 'UnitAveragePower', true);
                    case 'BFSK'
                        tx_sym = fskmod(bitsIn, M, 2000, 2, app.Fs);
                    case 'ASK'
                        % ASK -> On-Off Keying (1=amp, 0=0)
                        tx_sym = double(bitsIn); 
                end
                app.TxSym = tx_sym;
                
                % ---- Channel (AWGN) ----
                rx_sym = awgn(tx_sym, app.SnrSlider.Value, 'measured');
                
                % ---- Demodulation ----
                switch app.ModType
                    case 'BPSK', rx_bits_raw = pskdemod(rx_sym, M);
                    case 'QPSK', rx_bits_raw = pskdemod(rx_sym, M, pi/4, 'OutputType', 'bit');
                    case '16-QAM', rx_bits_raw = qamdemod(rx_sym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
                    case 'BFSK', rx_bits_raw = fskdemod(rx_sym, M, 2000, 2, app.Fs);
                    case 'ASK'
                        % ASK Demod -> Threshold at 0.5
                        rx_bits_raw = double(real(rx_sym) > 0.5);
                end
                if modPad>0, rx_bits_raw = rx_bits_raw(1:end-modPad); end
                app.RxBits_Ch = double(rx_bits_raw(:));
                
                % ---- Channel decoding ----
                if app.chkChannel.Value
                    n=7; k=4;
                    Lr = length(app.RxBits_Ch); rem7 = mod(Lr,7);
                    if rem7>0, app.RxBits_Ch = app.RxBits_Ch(1:end-rem7); end
                    if ~isempty(app.RxBits_Ch)
                        rx_words_mat = reshape(app.RxBits_Ch, n, []).';
                        decoded_words = decode(rx_words_mat, n, k, 'hamming/binary');
                        rx_bits_padded = reshape(decoded_words.', [], 1);
                        if app.HammingPadding>0
                            if length(rx_bits_padded) > app.HammingPadding
                                app.RxBits_Src = rx_bits_padded(1:end-app.HammingPadding);
                            else
                                app.RxBits_Src = [];
                            end
                        else
                            app.RxBits_Src = rx_bits_padded;
                        end
                    else
                        app.RxBits_Src = [];
                    end
                else
                    app.RxBits_Src = app.RxBits_Ch;
                end
                app.RxBits_Src = double(app.RxBits_Src(:));
                
                % ---- Source decoding (Huffman fallback robust) ----
                rx_sym_rec = []; use_pcm_fallback = false;
                if app.chkSource.Value
                    try
                        valid_len = length(app.TxBits_Src);
                        if length(app.RxBits_Src) >= valid_len, bits_to_decode = app.RxBits_Src(1:valid_len); else bits_to_decode = app.RxBits_Src; end
                        rx_sym_rec = huffmandeco(bits_to_decode, app.HuffDict);
                        if isempty(rx_sym_rec) || length(rx_sym_rec) ~= length(symbols), use_pcm_fallback = true; end
                    catch
                        use_pcm_fallback = true;
                    end
                else
                    use_pcm_fallback = true;
                end
                if use_pcm_fallback
                    L_bits = length(app.RxBits_Src);
                    nFullSymbols = floor(L_bits / qBit);
                    rx_sym_rec = [];
                    if nFullSymbols > 0
                        bits_to_dec = app.RxBits_Src(1:nFullSymbols*qBit);
                        rx_mat = reshape(bits_to_dec, qBit, []).';
                        rx_sym_rec = bi2de(rx_mat, 'left-msb');
                    end
                    needed = length(symbols) - length(rx_sym_rec);
                    if needed > 0
                        rnd_fill = randi([0 L-1], needed, 1);
                        rx_sym_rec = [rx_sym_rec; rnd_fill];
                    elseif needed < 0
                        rx_sym_rec = rx_sym_rec(1:length(symbols));
                    end
                end
                rx_sym_rec = double(rx_sym_rec(:));
                
                app.RxSymbols = rx_sym_rec; % Store for visualization
                
                rx_double = (double(rx_sym_rec) / (L-1)) * 2 - 1;
                N_orig = length(app.OriginalAudio); N_rx = length(rx_double);
                if N_rx < N_orig
                    needed = N_orig - N_rx;
                    noise_levels = randi([0 L-1], needed, 1);
                    noise_vals = (double(noise_levels) / (L-1)) * 2 - 1;
                    rx_double = [rx_double; noise_vals];
                elseif N_rx > N_orig
                    rx_double = rx_double(1:N_orig);
                end
                app.ProcessedAudio = rx_double;
                app.PlayRxBtn.Enable = 'on';
                if use_pcm_fallback && app.chkSource.Value
                    app.StatusLabel.Text = 'Status: Done (Noisy Fallback)';
                else
                    app.StatusLabel.Text = 'Status: Done';
                end
                
                % Plotting
                app.plotADC(x, x_int, raw_bin_strs, qBit);
                app.plotCodingTables();
                app.plotModulationWaves();
                app.plotDemodulationWaves(rx_sym);
                app.plotDecodingTables();
                app.plotOutput(rx_double);
            catch ME
                app.StatusLabel.Text = 'Error';
                uialert(app.Fig, ME.message, 'Sim Error');
            end
        end
        
        %% Plotting helpers
        function plotSource(app)
            delete(app.TabSource.Children);
            g = uigridlayout(app.TabSource, [1,1]);
            ax1 = uiaxes(g);
            plot(ax1, app.OriginalAudio); title(ax1,'Input Audio'); grid(ax1,'on');
        end
        
        function plotADC(app, analogSig, quantLevels, binStrs, bitsPerSample)
            delete(app.TabADC.Children);
            g = uigridlayout(app.TabADC, [1,1]);
            uit = uitable(g);
            N_table = min(200, length(analogSig));
            data = cell(N_table, 4);
            for i = 1:N_table
                data{i,1} = i;
                data{i,2} = analogSig(i);
                data{i,3} = quantLevels(i);
                if size(binStrs,1) >= i, data{i,4} = binStrs(i,:); else data{i,4} = repmat('0', 1, bitsPerSample); end
            end
            uit.Data = data;
            uit.ColumnName = {'Sample', 'Analog', 'Level', 'Binary'};
            uit.ColumnWidth = {60, 120, 80, '1x'};
        end
        
        function plotCodingTables(app)
            delete(app.TabCoding.Children);
            g = uigridlayout(app.TabCoding, [1,2]); 
            g.ColumnWidth = {'1x','1x'};
            
            % --- Panel 1: Source Coding (Split) ---
            p1 = uipanel(g,'Title','Source Coding');
            p1_g = uigridlayout(p1,[4,1]); 
            % FIXED HEADER HEIGHT (25px), FLEXIBLE TABLE HEIGHT (1x)
            p1_g.RowHeight = {25, '1x', 25, '1x'}; 
            p1_g.Padding = [5 5 5 5]; p1_g.RowSpacing = 5;
            
            % 1A: Reference Map
            l1 = uilabel(p1_g,'Text','1. Reference Dictionary (Symbol->Code)','FontWeight','bold');
            l1.Layout.Row = 1;
            
            tRef = uitable(p1_g,'RowStriping','on');
            tRef.Layout.Row = 2;
            
            if app.chkSource.Value && ~isempty(app.HuffDict)
                d = app.HuffDict;
                lens = cellfun(@length, d(:,2)); [~,idx] = sort(lens); d=d(idx,:);
                refData = cell(size(d,1),2);
                for i=1:size(d,1), refData{i,1}=d{i,1}; refData{i,2}=num2str(d{i,2}); end
                tRef.Data = refData; tRef.ColumnName = {'Symbol', 'Huffman Code'};
            else
                if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), qBit=4; else qBit=8; end
                L = 2^qBit;
                refData = cell(L,2);
                bin = dec2bin(0:L-1, qBit);
                for i=1:L, refData{i,1}=i-1; refData{i,2}=bin(i,:); end
                tRef.Data = refData; tRef.ColumnName = {'Symbol', 'PCM Code'};
            end
            % 1B: Real-Time Stream
            l2 = uilabel(p1_g,'Text','2. Real-Time Stream (Tx)','FontWeight','bold');
            l2.Layout.Row = 3;
            
            tStream = uitable(p1_g,'RowStriping','on');
            tStream.Layout.Row = 4;
            
            if ~isempty(app.TxSymbols)
                numShow = min(50, length(app.TxSymbols));
                sData = cell(numShow,3);
                
                if app.chkSource.Value && ~isempty(app.HuffDict)
                    hDict = app.HuffDict; hSyms = [hDict{:,1}];
                    for i=1:numShow
                        val = app.TxSymbols(i);
                        sData{i,1} = i; sData{i,2} = val;
                        idx = find(hSyms == val, 1);
                        if ~isempty(idx), sData{i,3} = num2str(hDict{idx,2}); else sData{i,3}='ERR'; end
                    end
                else
                    if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), qBit=4; else qBit=8; end
                    bin = dec2bin(app.TxSymbols(1:numShow), qBit);
                    for i=1:numShow, sData{i,1}=i; sData{i,2}=app.TxSymbols(i); sData{i,3}=bin(i,:); end
                end
                tStream.Data = sData; tStream.ColumnName = {'Sample', 'Symbol', 'Tx Code'};
                tStream.ColumnWidth = {60, 60, '1x'};
            end
            
            % --- Panel 2: Channel Coding (Split) ---
            p2 = uipanel(g,'Title','Channel Coding (Hamming)');
            p2_grid = uigridlayout(p2,[4,1]); 
            p2_grid.RowHeight = {25, '1x', 25, '1x'};
            p2_grid.Padding = [5 5 5 5]; p2_grid.RowSpacing = 5;
            
            % 2A: Reference Map
            l3 = uilabel(p2_grid,'Text','1. Reference Mapping (Static)','FontWeight','bold');
            l3.Layout.Row = 1;
            
            tRefCh = uitable(p2_grid,'RowStriping','on');
            tRefCh.Layout.Row = 2;
            
            if app.chkChannel.Value
                in_ref = de2bi((0:15)', 4, 'left-msb');
                out_ref = encode(in_ref, 7, 4, 'hamming/binary');
                chData = cell(16,2);
                for i=1:16, chData{i,1}=num2str(in_ref(i,:)); chData{i,2}=num2str(out_ref(i,:)); end
                tRefCh.Data = chData; tRefCh.ColumnName = {'Data(4)', 'Code(7)'};
                tRefCh.ColumnWidth = {'1x', '1x'};
            else
                uilabel(p2_grid,'Text','Disabled'); tRefCh.Visible='off';
            end
            % 2B: Real-Time Stream
            l4 = uilabel(p2_grid,'Text','2. Real-Time Stream (Tx)','FontWeight','bold');
            l4.Layout.Row = 3;
            
            tRealCh = uitable(p2_grid,'RowStriping','on');
            tRealCh.Layout.Row = 4;
            
            if app.chkChannel.Value && ~isempty(app.TxBits_Src)
                nBlocks=50;
                if length(app.TxBits_Src) >= 4*nBlocks
                    in_mat = reshape(app.TxBits_Src(1:4*nBlocks), 4, []).';
                    out_mat = encode(in_mat, 7, 4, 'hamming/binary');
                    hData = cell(nBlocks,2);
                    for i=1:nBlocks, hData{i,1}=num2str(in_mat(i,:)); hData{i,2}=num2str(out_mat(i,:)); end
                    tRealCh.Data = hData; tRealCh.ColumnName = {'Stream In', 'Stream Out'};
                    tRealCh.ColumnWidth = {'1x', '1x'};
                else
                    tRealCh.Visible='off';
                end
            else
                tRealCh.Visible='off';
            end
        end
        
        function plotModulationWaves(app)
            delete(app.TabModulation.Children);
            g = uigridlayout(app.TabModulation, [3,1]);
            if isempty(app.TxBits_Ch), return; end
            bitsToShow = app.TxBits_Ch(1:min(20,end)); bitsToShow = bitsToShow(:);
            Rb = 100; Fc = 400; Fs_plot = 10000; samplesPerBit = Fs_plot/Rb;
            ax1 = uiaxes(g);
            bit_wave_col = kron(bitsToShow, ones(samplesPerBit,1));
            total_samples = length(bit_wave_col); t_vec = (0:total_samples-1)*(1/Fs_plot); t_row = t_vec(:).'; b_row = bit_wave_col(:).';
            plot(ax1, t_row, b_row, 'k', 'LineWidth', 2); title(ax1,'Tx Bits'); ylim(ax1,[-0.2 1.2]); grid(ax1,'on');
            ax2 = uiaxes(g); carrier_wave = cos(2*pi*Fc*t_row); plot(ax2,t_row,carrier_wave,'Color',[0.6 0.6 0.6]); title(ax2,'Carrier'); ylim(ax2,[-1.2 1.2]); grid(ax2,'on');
            ax3 = uiaxes(g);
            mod_wave = zeros(size(t_row));
            switch app.ModType
                case 'BPSK'
                    mod_wave = (2*b_row - 1) .* carrier_wave;
                case 'BFSK'
                    f1=Fc; f2=Fc*2;
                    for i=1:length(bitsToShow)
                        idx_start = (i-1)*samplesPerBit+1; idx_end = i*samplesPerBit;
                        if idx_end>length(t_row), break; end
                        idx = idx_start:idx_end; t_seg = t_row(idx);
                        if bitsToShow(i)==0, mod_wave(idx)=cos(2*pi*f1*t_seg); else mod_wave(idx)=cos(2*pi*f2*t_seg); end
                    end
                case 'ASK'
                    % OOK (Bit 1 = Carrier, Bit 0 = Flat)
                    mod_wave = b_row .* carrier_wave;
                otherwise
                    mod_wave = (2*b_row - 1).*carrier_wave;
            end
            plot(ax3, t_row, mod_wave, 'b'); title(ax3, ['Modulated: ' app.ModType]); grid(ax3,'on');
            linkaxes([ax1 ax2 ax3], 'x'); xlim(ax1, [0 t_row(end)]);
        end
        
        function plotDemodulationWaves(app, rx_sym)
            delete(app.TabDemod.Children);
            g = uigridlayout(app.TabDemod, [2,1]); g.RowHeight = {'2x','1x'};
            if isempty(app.RxBits_Ch), return; end
            
            % --- Top part: Waveforms ---
            p1 = uipanel(g,'Title','Physical Layer Vis'); g1 = uigridlayout(p1,[3,1]); g1.Padding=[0 0 0 0]; g1.RowSpacing=2;
            bitsToShow = app.RxBits_Ch(1:min(20,end)); bitsToShow = bitsToShow(:);
            Rb=100; Fc=400; Fs_plot=10000; samplesPerBit=Fs_plot/Rb;
            total_bits = length(bitsToShow); total_samples = total_bits * samplesPerBit;
            t_row = (0:total_samples-1)*(1/Fs_plot);
            mod_wave_tx = zeros(size(t_row));
            switch app.ModType
                case 'BPSK'
                    bit_wave_ref = kron(bitsToShow, ones(samplesPerBit,1));
                    mod_wave_tx = (2*bit_wave_ref(:).' - 1) .* cos(2*pi*Fc*t_row);
                case 'BFSK'
                    for i=1:length(bitsToShow)
                        idx_start = (i-1)*samplesPerBit+1; idx_end = i*samplesPerBit;
                        if idx_end>length(t_row), idx_end=length(t_row); end
                        if bitsToShow(i)==0, mod_wave_tx(idx_start:idx_end)=cos(2*pi*Fc*t_row(idx_start:idx_end));
                        else mod_wave_tx(idx_start:idx_end)=cos(2*pi*(Fc*2)*t_row(idx_start:idx_end)); end
                    end
                case 'ASK'
                    % OOK
                    bit_wave_ref = kron(bitsToShow, ones(samplesPerBit,1));
                    mod_wave_tx = bit_wave_ref(:).' .* cos(2*pi*Fc*t_row);
                otherwise
                    bit_wave_ref = kron(bitsToShow, ones(samplesPerBit,1));
                    mod_wave_tx = (2*bit_wave_ref(:).' - 1) .* cos(2*pi*Fc*t_row);
            end
            noisy_wave = awgn(mod_wave_tx, app.SnrSlider.Value, 'measured');
            axTx = uiaxes(g1); plot(axTx, t_row, mod_wave_tx, 'Color', [0.5 0.5 0.5]); title(axTx,'Tx Wave');
            axRx = uiaxes(g1); plot(axRx, t_row, noisy_wave, 'b'); title(axRx,'Rx Wave');
            axBits = uiaxes(g1); stairs(axBits, bitsToShow, 'r', 'LineWidth', 2); title(axBits,'Recovered Bits'); ylim(axBits,[-0.2 1.2]);
            linkaxes([axTx axRx],'x');
            
            % --- Bottom part: Constellations ---
            % Create a 1x2 grid for the two constellation plots
            g2 = uigridlayout(g, [1,2]);
            
            % Tx Constellation
            axTxConst = uiaxes(g2);
            if ~isempty(app.TxSym)
                % Plot a subset of points for performance
                numPoints = min(500, length(app.TxSym));
                scatter(axTxConst, real(app.TxSym(1:numPoints)), imag(app.TxSym(1:numPoints)), 'b.');
            end
            title(axTxConst, 'Tx Constellation'); 
            axis(axTxConst,'equal');
            grid(axTxConst, 'on');
            
            % Rx Constellation
            axRxConst = uiaxes(g2);
            numPoints = min(500, length(rx_sym));
            scatter(axRxConst, real(rx_sym(1:numPoints)), imag(rx_sym(1:numPoints)), 'r.');
            title(axRxConst, 'Rx Constellation'); 
            axis(axRxConst,'equal');
            grid(axRxConst, 'on');
        end
        
        function plotDecodingTables(app)
            delete(app.TabDecoding.Children);
            g = uigridlayout(app.TabDecoding, [1,2]); 
            g.ColumnWidth = {'1x','1x'};
            
            % --- Panel 1: Channel Decoding (Split) ---
            p1 = uipanel(g,'Title','Channel Decoding (Hamming)');
            p1_grid = uigridlayout(p1,[4,1]);
            p1_grid.RowHeight = {25, '1x', 25, '1x'};
            p1_grid.Padding = [5 5 5 5]; p1_grid.RowSpacing = 5;
            
            % 1A: Reference
            l1 = uilabel(p1_grid, 'Text', '1. Reference (Codeword -> Data)', 'FontWeight', 'bold');
            l1.Layout.Row = 1;
            
            tRefCh = uitable(p1_grid,'RowStriping','on');
            tRefCh.Layout.Row = 2;
            
            if app.chkChannel.Value
                in_ref = de2bi((0:15)', 4, 'left-msb');
                out_ref = encode(in_ref, 7, 4, 'hamming/binary');
                refData = cell(16,2);
                for i=1:16
                    refData{i,1} = num2str(out_ref(i,:)); 
                    refData{i,2} = num2str(in_ref(i,:)); 
                end
                tRefCh.Data = refData; tRefCh.ColumnName = {'Codeword(7)', 'Data(4)'};
                tRefCh.ColumnWidth = {'1x','1x'};
            else
                tRefCh.Visible='off'; uilabel(p1_grid,'Text','Disabled');
            end
            % 1B: Stream
            l2 = uilabel(p1_grid, 'Text', '2. Real-Time Stream (Rx)', 'FontWeight', 'bold');
            l2.Layout.Row = 3;
            
            tRealCh = uitable(p1_grid,'RowStriping','on');
            tRealCh.Layout.Row = 4;
            
            if app.chkChannel.Value && ~isempty(app.RxBits_Ch)
                nBlocks=50;
                if length(app.RxBits_Ch) >= 7*nBlocks
                    rx_mat = reshape(app.RxBits_Ch(1:7*nBlocks), 7, []).';
                    dec_mat = decode(rx_mat, 7, 4, 'hamming/binary');
                    cData = cell(nBlocks,2);
                    for i=1:nBlocks, cData{i,1}=num2str(rx_mat(i,:)); cData{i,2}=num2str(dec_mat(i,:)); end
                    tRealCh.Data = cData; tRealCh.ColumnName = {'Rx(7)', 'Decoded(4)'};
                    tRealCh.ColumnWidth = {'1x','1x'};
                else
                    tRealCh.Visible='off';
                end
            else
                 tRealCh.Visible='off';
            end
            
            % --- Panel 2: Source Decoding (Split) ---
            p2 = uipanel(g,'Title','Source Decoding');
            p2_grid = uigridlayout(p2,[4,1]);
            p2_grid.RowHeight = {25, '1x', 25, '1x'};
            p2_grid.Padding = [5 5 5 5]; p2_grid.RowSpacing = 5;
            
            % 2A: Reference
            l3 = uilabel(p2_grid,'Text','1. Reference (Code->Symbol)','FontWeight','bold');
            l3.Layout.Row = 1;
            
            tRefSrc = uitable(p2_grid,'RowStriping','on');
            tRefSrc.Layout.Row = 2;
            
            if app.chkSource.Value && ~isempty(app.HuffDict)
                d = app.HuffDict;
                lens = cellfun(@length, d(:,2)); [~,idx] = sort(lens); d=d(idx,:);
                refData = cell(size(d,1),2);
                for i=1:size(d,1), refData{i,1}=num2str(d{i,2}); refData{i,2}=d{i,1}; end
                tRefSrc.Data = refData; tRefSrc.ColumnName = {'Code', 'Symbol'};
                tRefSrc.ColumnWidth = {'1x', 60};
            else
                if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), qBit=4; else qBit=8; end
                L = 2^qBit;
                refData = cell(L,2);
                bin = dec2bin(0:L-1, qBit);
                for i=1:L, refData{i,1}=bin(i,:); refData{i,2}=i-1; end
                tRefSrc.Data = refData; tRefSrc.ColumnName = {'PCM Code', 'Symbol'};
                tRefSrc.ColumnWidth = {'1x', 60};
            end
            
            % 2B: Stream
            l4 = uilabel(p2_grid,'Text','2. Real-Time Stream (Rx)','FontWeight','bold');
            l4.Layout.Row = 3;
            
            tRealSrc = uitable(p2_grid,'RowStriping','on');
            tRealSrc.Layout.Row = 4;
            
            if ~isempty(app.RxSymbols)
                numShow = min(50, length(app.RxSymbols));
                if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), L=16; else L=256; end
                sData = cell(numShow,3);
                
                if app.chkSource.Value && ~isempty(app.HuffDict)
                    hDict = app.HuffDict; hSyms = [hDict{:,1}];
                    for i=1:numShow
                        val = app.RxSymbols(i);
                        sData{i,1} = val;
                        idx = find(hSyms == val, 1);
                        if ~isempty(idx), sData{i,2} = num2str(hDict{idx,2}); else sData{i,2}='?'; end
                        sData{i,3} = sprintf('%.3f', (double(val)/(L-1))*2 - 1);
                    end
                else
                    if L==16, qBit=4; else qBit=8; end
                    bin = dec2bin(app.RxSymbols(1:numShow), qBit);
                    for i=1:numShow
                        val = app.RxSymbols(i);
                        sData{i,1} = val;
                        sData{i,2} = bin(i,:);
                        sData{i,3} = sprintf('%.3f', (double(val)/(L-1))*2 - 1);
                    end
                end
                tRealSrc.Data = sData; tRealSrc.ColumnName = {'Symbol', 'Code', 'Analog'};
                tRealSrc.ColumnWidth = {60, '1x', 80};
            else
                tRealSrc.Visible='off';
            end
        end
        
        function plotOutput(app, rx)
            delete(app.TabOutput.Children);
            g = uigridlayout(app.TabOutput, [3,1]); g.RowHeight = {'2x', '1x', '1x'};
            
            % 1. Audio Overlay
            ax1 = uiaxes(g);
            t = (0:length(app.OriginalAudio)-1)/app.Fs;
            plot(ax1, t, app.OriginalAudio, 'b'); hold(ax1,'on');
            plot(ax1, t, rx, 'r--'); legend(ax1, 'Tx', 'Rx'); title(ax1,'Audio Overlay'); grid(ax1,'on');
            
            % 2. Info panel
            axInfo = uiaxes(g); axInfo.Visible = 'off';
            L_ber = min(length(app.TxBits_Ch), length(app.RxBits_Ch));
            if L_ber>0
                [nErr, ber] = biterr(app.TxBits_Ch(1:L_ber), app.RxBits_Ch(1:L_ber));
            else
                nErr = 0; ber = NaN;
            end
            info = {
                sprintf('Coding: Src=%d, Ch=%d', app.chkSource.Value, app.chkChannel.Value),
                sprintf('Modulation: %s', app.ModType),
                sprintf('Tx Bits: %d', length(app.TxBits_Ch)),
                sprintf('Bit Errors: %d', nErr),
                sprintf('BER (instant): %g', ber),
                '',
                sprintf('Status: %s', app.StatusLabel.Text)
            };
            text(axInfo, 0.1, 0.5, info, 'FontSize', 12, 'FontName', 'Courier');
            
            % 3. BER vs SNR axes
            axBer = uiaxes(g); grid(axBer,'on');
            title(axBer, ['BER vs SNR (' app.ModType ')']);
            xlabel(axBer, 'SNR (dB)'); ylabel(axBer, 'BER'); set(axBer,'YScale','log');
            
            % A. Theoretical Curve
            SNRdB_th = -5:1:25;
            ber_th = [];
            switch app.ModType
                case 'BPSK', ber_th = berawgn(SNRdB_th, 'psk', 2, 'nondiff');
                case 'QPSK', ber_th = berawgn(SNRdB_th, 'psk', 4, 'nondiff');
                case '16-QAM', ber_th = berawgn(SNRdB_th, 'qam', 16);
                case 'BFSK', ber_th = berawgn(SNRdB_th, 'fsk', 2, 'coherent');
                case 'ASK'
                    % OOK approx: unipolar is 3dB worse than BPSK
                    EbN0 = 10.^(SNRdB_th./10);
                    ber_th = 0.5 * erfc(sqrt(EbN0 ./ 4));
            end
            hold(axBer, 'on'); 
            if ~isempty(ber_th)
                plot(axBer, SNRdB_th, ber_th, 'k-', 'LineWidth',1.5, 'DisplayName','Theory');
            end
            
            % B. Simulated Curve
            SNRdB_sim = -5:2:25;
            if isempty(app.TxBits_Ch)
                tmpTxBits = app.generateTxBitsForBer();
            else
                tmpTxBits = app.TxBits_Ch;
            end
            
            if ~isempty(tmpTxBits)
                simBer = app.simulateBerVsSnr(tmpTxBits, SNRdB_sim);
                plot(axBer, SNRdB_sim, simBer, 'ro-', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
            else
                text(axBer, 0.2, 0.5, 'No Tx bits', 'Units', 'normalized');
            end
            
            % C. CURRENT SNR Marker (The Green Dot)
            currentSnr = app.SnrSlider.Value;
            currentBerTh = 1e-6; 
            switch app.ModType
                case 'BPSK', currentBerTh = berawgn(currentSnr, 'psk', 2, 'nondiff');
                case 'QPSK', currentBerTh = berawgn(currentSnr, 'psk', 4, 'nondiff');
                case '16-QAM', currentBerTh = berawgn(currentSnr, 'qam', 16);
                case 'BFSK', currentBerTh = berawgn(currentSnr, 'fsk', 2, 'coherent');
                case 'ASK'
                    ebn0_curr = 10^(currentSnr/10);
                    currentBerTh = 0.5 * erfc(sqrt(ebn0_curr / 4));
            end
            plot(axBer, currentSnr, currentBerTh, 'h', ...
                'MarkerSize', 12, 'MarkerFaceColor', 'g', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
                'DisplayName', sprintf('Current SNR (%.1f dB)', currentSnr));
                
            legend(axBer, 'Location', 'southwest');
            ylim(axBer, [1e-6 1]);
        end
        
        %% BER helper: regenerate TxBits_Ch from current audio/settings (non-destructive)
        function txBits = generateTxBitsForBer(app)
            txBits = [];
            try
                max_amp = max(abs(app.OriginalAudio));
                if max_amp > 1e-6, x = app.OriginalAudio / max_amp; else x = app.OriginalAudio; end
                if strcmp(app.QuantizationDropDown.Value, '4-bit (16 levels)'), qBit = 4; else qBit = 8; end
                L = 2^qBit;
                x_int = floor(((x + 1)/2)*(L-1));
                x_int(x_int<0)=0; x_int(x_int> L-1)=L-1;
                symbols = double(x_int);
                
                if app.chkSource.Value
                    tabS = tabulate(symbols);
                    unique_sym = tabS(tabS(:,2)>0,1);
                    probs = tabS(tabS(:,2)>0,3)/100;
                    hdict = huffmandict(unique_sym, probs);
                    bits_src = huffmanenco(symbols, hdict);
                else
                    bits_src = reshape(de2bi(symbols, qBit, 'left-msb').', [], 1);
                end
                bits_src = double(bits_src(:));
                
                if app.chkChannel.Value
                    n=7; k=4;
                    remBits = mod(length(bits_src), k); padding = 0;
                    if remBits>0, padding = k-remBits; bits_in = [bits_src; zeros(padding,1)]; else bits_in = bits_src; end
                    msg_words = reshape(bits_in, k, []).';
                    coded_words = encode(msg_words, n, k, 'hamming/binary');
                    txBits = reshape(coded_words.', [], 1);
                else
                    txBits = bits_src;
                end
                txBits = double(txBits(:));
            catch
                txBits = [];
            end
        end
        
        %% Simulate BER vs SNR for given bitstream using current modulation config
        function simBer = simulateBerVsSnr(app, txBits, SNRdB_vec)
            if isempty(txBits), simBer = []; return; end
            M = 2;
            switch app.ModType
                case 'BPSK', M=2;
                case 'QPSK', M=4;
                case '16-QAM', M=16;
                case 'BFSK', M=2;
                case 'ASK', M=2; % OOK
            end
            bps = log2(M);
            bitsIn = txBits;
            modPad = 0;
            if mod(length(bitsIn), bps) ~= 0
                modPad = bps - mod(length(bitsIn), bps);
                bitsIn = [bitsIn; zeros(modPad, 1)];
            end
            
            % Modulate once (OOK handled here)
            switch app.ModType
                case 'BPSK'
                    txSym = pskmod(bitsIn, M);
                case 'QPSK'
                    txSym = pskmod(bitsIn, M, pi/4, 'InputType', 'bit');
                case '16-QAM'
                    txSym = qammod(bitsIn, M, 'InputType', 'bit', 'UnitAveragePower', true);
                case 'BFSK'
                     if bps>1
                        bitsMat = reshape(bitsIn, bps, []).';
                        intSyms = bi2de(bitsMat, 'left-msb');
                    else
                        intSyms = bitsIn;
                    end
                    txSym = fskmod(intSyms, M, 2000, 2, app.Fs);
                case 'ASK'
                    % OOK
                    txSym = double(bitsIn);
                otherwise
                    txSym = pskmod(bitsIn, M);
            end
            
            simBer = zeros(size(SNRdB_vec));
            for ii = 1:length(SNRdB_vec)
                SNRdB = SNRdB_vec(ii);
                rxSym = awgn(txSym, SNRdB, 'measured');
                switch app.ModType
                    case 'BPSK'
                        rx_bits_raw = pskdemod(rxSym, M);
                    case 'QPSK'
                        rx_bits_raw = pskdemod(rxSym, M, pi/4, 'OutputType', 'bit');
                    case '16-QAM'
                        rx_bits_raw = qamdemod(rxSym, M, 'OutputType', 'bit', 'UnitAveragePower', true);
                    case 'BFSK'
                        rx_int = fskdemod(rxSym, M, 2000, 2, app.Fs);
                         if bps>1
                            rx_bits_mat = de2bi(rx_int, bps, 'left-msb');
                            rx_bits_raw = reshape(rx_bits_mat.', [], 1);
                        else
                            rx_bits_raw = rx_int(:);
                        end
                    case 'ASK'
                        % OOK Demod
                        rx_bits_raw = double(real(rxSym) > 0.5);
                    otherwise
                        rx_bits_raw = pskdemod(rxSym, M);
                end
                
                if modPad>0 && length(rx_bits_raw) >= modPad
                    rx_bits_raw = rx_bits_raw(1:end-modPad);
                end
                rx_bits_raw = double(rx_bits_raw(:));
                Lc = min(length(txBits), length(rx_bits_raw));
                if Lc>0
                    [~, ber_val] = biterr(txBits(1:Lc), rx_bits_raw(1:Lc));
                    simBer(ii) = ber_val;
                else
                    simBer(ii) = NaN;
                end
            end
        end
        
        %% Play functions
        function playOriginal(app, ~, ~)
            if ~isempty(app.OriginalAudio), soundsc(app.OriginalAudio, app.Fs); end
        end
        
        function playProcessed(app, ~, ~)
            if ~isempty(app.ProcessedAudio), soundsc(app.ProcessedAudio, app.Fs); end
        end
    end
end