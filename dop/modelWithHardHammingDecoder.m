function [PeBit, PED, T] = modelWithHardHammingDecoder(k, gX, codes, sigma, N, modulSignalBook, matrixHammingDec)

r = length(gX) - 1;
n = k + r + 5;
K = 2^k;

Ncur = 0;
Nt = 0;

PeBit = 0;
PED = 0;

%     ��������
indexCode = randi([1 K],1,1);%��������� ������� ���. �����

while Ncur < N
    mS = modulSignalBook(indexCode, :);
    
    %     ����
    mR = mS + sigma * randn(1, n);
    
    %     BPSK ^-1
    mX_ = mR < 0;
    
    %     Hamming Decoder
    extended_mX_ = [0, 0, mX_]; % ��������� ��� ���� � ������
    extended_mH_ = hardHammingDecoder(matrixHammingDec, extended_mX_);
    mH_ = extended_mH_(3:1:end); % ������� ��� ���� � ������
    
    %     CRC-r ^-1    
    flagCRC = sum(modGx(mH_, gX));
    flagSum = sum(xor(mH_, codes(indexCode, :)));
    
    if flagCRC == 0 % ������ ��� ��� �� ����������
        PED = PED + (flagSum > 0 & flagCRC == 0);
        indexCode = randi([1 K],1,1);%��������� ������� ���. �����
        Ncur = Ncur + 1;
    end
    
    PeBit = PeBit + flagSum;
    Nt = Nt + 1;
end

PeBit = PeBit / Nt / n;
PED = PED / Nt;
T = (k * N) / (n * Nt);