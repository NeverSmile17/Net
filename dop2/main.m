clc
clear
close

E = 1;
N = 100000;
k = 4;
r = 3;
n = k + r;

[informBook, modulSignalBook] = modulatedCodeBook(k);

SNRdB = -20 : 0;
PeBits_1 = zeros (1, length(SNRdB));
PeBits_2 = zeros (1, length(SNRdB));

for i = 1 : length(SNRdB)
    
    disp(SNRdB(i));
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(E / (2*SNR));
    
    PeBits_1(1, i) = model_1(informBook, modulSignalBook, k, n, sigma, N);
    PeBits_2(1, i) = model_2(informBook, modulSignalBook, k, n, sigma, N);
end

figure();
axis('square');
semilogy(SNRdB, PeBits_1, 'b.-', SNRdB, PeBits_2, 'r-');
xlabel('SNRdB'); 
ylabel('PeBit');
legend ({'Вероятность ошибки на бит при жёстком декодировании', ...
         'Вероятность ошибки на бит при мягком декодировании'}, ...
'Location','southwest')

%Кодер Хэмминга

function mH = HammingEncoder(matrixHammingEnc, mX)	
    mH = mod(mX*matrixHammingEnc, 2);
end

%Мягкое декодирование Хэмминга

function mX = softHammingDecoder(informBook, modulSignalBook, mR)
    [~, I] = min(pdist2(modulSignalBook(1:end, :), mR));
    mX = informBook(I, :);
end

%Жесткое декодирование Хэмминга

function mX = hardHammingDecoder(matrixHammingDec, mH)
    S = mod(mH * matrixHammingDec, 2);
    if sum(S) == 0
        mX = mH(1:1:end - 3);
    else
        indexError = bi2de(S);
        mH(indexError) = xor(mH(indexError), 1);
        mX = mH(1:1:end - 3);
    end
end