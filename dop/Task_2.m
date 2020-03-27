clc
clear all
close all

%                   1 2 3 4 5
matrixHammingDec = [1 0 0 0 0; % 1
                    0 1 0 0 0; % 2
                    1 1 0 0 0; % 3
                    0 0 1 0 0; % 4
                    1 0 1 0 0; % 5
                    0 1 1 0 0; % 6
                    1 1 1 0 0; % 7
                    0 0 0 1 0; % 8
                    1 0 0 1 0; % 9
                    0 1 0 1 0; % 10
                    1 1 0 1 0; % 11
                    0 0 1 1 0; % 12
                    1 0 1 1 0; % 13
                    0 1 1 1 0; % 14
                    1 1 1 1 0; % 15
                    0 0 0 0 1; % 16
                    1 0 0 0 1; % 17
                    0 1 0 0 1; % 18
                    1 1 0 0 1; % 19
                    0 0 1 0 1; % 20
                    1 0 1 0 1; % 21
                    0 1 1 0 1; % 22
                    1 1 1 0 1; % 23
                    0 0 0 1 1; % 24
                    1 0 0 1 1; % 25
                    0 1 0 1 1; % 26
                    % % % % %
                    1 0 0 0 0; % 27
                    0 1 0 0 0; % 28
                    0 0 1 0 0; % 29
                    0 0 0 1 0; % 30
                    0 0 0 0 1];% 31

E = 1;
k = 8;
N = 100000; % ���������� ��������� ��� �������������
               
%gX = x^16+x^13+x^12+x^11+x^10+x^8+x^6+x^5+x^2+1
gX = [1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1];
r = length(gX) - 1;
n = r + k;

codes = codeBook(k, gX); % ������� �����
modulSignalBook = modulatedCodeBook(codes); % ����������������� ������� ����� � ���������������� �����

d_min = min(sum(codes(2:end, :),2));
A = A_func(codes);    ����� ������� ���� ������� ����   

SNRdB = -8 : 2 : 16;
PeBits_hard = zeros (1, length(SNRdB));
PEDdB_hard = zeros (1, length(SNRdB));
TdB_hard = zeros (1, length(SNRdB));

PeBits_soft = zeros (1, length(SNRdB));
PEDdB_soft = zeros (1, length(SNRdB));
TdB_soft = zeros (1, length(SNRdB));

parfor i = 1 : length(SNRdB)
    tic
    
    disp(SNRdB(i));
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(E / (2*SNR));
    
    [PeBit_hard, PED_hard, T_hard] = modelWithHardHammingDecoder(k, gX,...
        codes, sigma, N, modulSignalBook, matrixHammingDec);
    
    PeBits_hard(1, i) = PeBit_hard;
    PEDdB_hard(1, i) = PED_hard;
    TdB_hard(1, i) = T_hard;
    
    [PeBit_soft, PED_soft, T_soft] = modelWithSoftHammingDecoder(k, gX,...
        codes, sigma, N, modulSignalBook);
    
    PeBits_soft(1, i) = PeBit_soft;
    PEDdB_soft(1, i) = PED_soft;
    TdB_soft(1, i) = T_soft;
    
    toc
end

%������ ����������� ������ ������������� CRC:
%������������� PED
SNRtheor = 10.^(SNRdB/10);
PeBitstheor = qfunc(sqrt(2*SNRtheor));
PEDsExact = zeros (1, length(SNRdB));

for i = 1 : length(SNRdB)
    for j = d_min : n
        PEDsExact(1, i) = PEDsExact(1, i) + A(j + 1) * ...
            PeBitstheor(i)^j * (1 - PeBitstheor(i))^(n - j);
    end
end

%�������� ���������� ����������� � ���������� ����
NoNoiseT = ones(1, length(SNRdB)).*k ./ (n + 5);

figure();
axis('square');
semilogy(SNRdB, PeBits_hard, 'b.-', SNRdB, PeBits_soft, 'cx-', SNRdB, PeBitstheor, 'r-');
xlabel('SNRdB'); 
ylabel('PeBit');
legend ({'�������� ����������� ������ �� ��� � �������������� ��������������� ����������� � ������ ���������', ...
    '�������� ����������� ������ �� ��� � �������������� ��������������� ����������� � ������ ���������', ...
    '�������� ����������� ������ �� ��� ��� ������������� ��������������� �����������'}, ...
'Location','southwest')

figure();
axis('square');
semilogy(SNRdB, PEDdB_hard, 'b.-', SNRdB, PEDsExact, 'r-');
xlabel('SNRdB'); 
ylabel('PED');
legend({'�������� ����������� ������ �������� � �������������� ��������������� ����������� � ������ ���������', ...
    '�������� ����������� ������ �������� ��� ������������� ��������������� �����������'},...
'Location','east')

figure();
axis('square');
hold on
semilogy(SNRdB, TdB_hard, 'b.-', SNRdB, NoNoiseT, 'r-');
xlabel('SNRdB'); 
ylabel('T');
legend ({'���������� ����������� ������� � �������������� ��������������� ����������� � ������ ��������� ��� ����', ...
    '���������� ����������� ������� � �������������� ��������������� ����������� ��� ����'}, ...
'Location','southwest')