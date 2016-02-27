classdef DuNucsTime
    properties (Constant)
        %A matrix 
        %Column 1: Index
        %Column 2: # Nuclei
        %Column 3: # E-cells at given time/nuclei
        findE = [1 2 0
                 2 2 0
                 3 2 0
                 4 3 0
                 5 4 0
                 6 4 0
                 7 4 0
                 8 4 0
                 9 4 0
                 10 4 0
                 11 4 0
                 12 4 0
                 13 6 0
                 14 6 0
                 15 7 1
                 16 7 1
                 17 7 1
                 18 8 1
                 19 8 1
                 20 8 1
                 21 8 1
                 22 8 1
                 23 11 1
                 24 12 1
                 25 12 1
                 26 13 1
                 27 14 2
                 28 14 2
                 29 14 2
                 30 15 2
                 31 15 2
                 32 15 2
                 33 15 2
                 34 16 2
                 35 24 2
                 36 24 2
                 37 24 2
                 38 24 2
                 39 26 2
                 40 26 2
                 41 26 2
                 42 26 2
                 43 26 2
                 44 26 2
                 45 28 2
                 46 28 2
                 47 28 2
                 48 28 2
                 49 28 2
                 50 31 2
                 51 42 2
                 52 45 4
                 53 46 4
                 54 46 4
                 55 46 4
                 56 47 4
                 57 48 4
                 58 51 4
                 59 51 4
                 60 51 4
                 61 51 4
                 62 51 4
                 63 51 4
                 64 51 4
                 65 53 4
                 66 53 4
                 67 54 4
                 68 55 4
                 69 57 4
                 70 76 4
                 71 85 4
                 72 86 4
                 73 87 4
                 74 87 4
                 75 88 4
                 76 89 4
                 77 93 5
                 78 97 6
                 79 98 6
                 80 98 6
                 81 99 7
                 82 100 8
                 83 100 8
                 84 100 8
                 85 100 8
                 86 101 8
                 87 102 8
                 88 102 8
                 89 102 8
                 90 104 8
                 91 109 8
                 92 115 8
                 93 127 8
                 94 137 8
                 95 149 8
                 96 153 8
                 97 160 8
                 98 167 8
                 99 172 8
                 100 176 8
                 101 177 8
                 102 177 8
                 103 178 8
                 104 182 8
                 105 184 8
                 106 186 8
                 107 186 8
                 108 187 8
                 109 188 8
                 110 188 8
                 111 190 8
                 112 190 8
                 113 190 8
                 114 190 8
                 115 190 8
                 116 192 8
                 117 193 8
                 118 193 8
                 119 196 8
                 120 200 8
                 121 201 8
                 122 205 8
                 123 212 8
                 124 227 8
                 125 240 8
                 126 252 9
                 127 273 10
                 128 282 11
                 129 290 11
                 130 300 11
                 131 313 12
                 132 321 13
                 133 329 13
                 134 333 13
                 135 339 13
                 136 347 14
                 137 352 15
                 138 354 16
                 139 356 16
                 140 356 16
                 141 356 16
                 142 359 16
                 143 361 16
                 144 361 16
                 145 361 16];
        
        
        
        %Time in min in first column, nuclei counts in second column   
        divTimes = [1 1
                    2 1
                    3 1
                    4 1
                    5 1
                    6 1
                    7 1
                    8 1
                    9 1
                    10 1
                    11 1
                    12 1
                    13 1
                    14 1
                    15 1
                    16 1
                    17 1
                    18 1
                    19 1
                    20 1
                    21 1
                    22 1
                    23 1
                    24 1
                    25 1
                    26 1
                    27 1
                    28 1
                    29 1
                    30 1
                    31 1
                    32 1
                    33 1
                    34 1
                    35 2
                    36 2
                    37 2
                    38 2
                    39 2
                    40 2
                    41 2
                    42 2
                    43 2
                    44 2
                    45 3
                    46 3
                    47 3
                    48 4
                    49 4
                    50 5
                    51 5
                    52 7
                    53 7
                    54 6
                    55 8
                    56 8
                    57 7
                    58 7
                    59 8
                    60 8
                    61 8
                    62 8
                    63 8
                    64 8
                    65 8
                    66 10
                    67 12
                    68 12
                    69 12
                    70 12
                    71 12
                    72 13
                    73 14
                    74 14
                    75 14
                    76 14
                    77 15
                    78 15
                    79 15
                    80 15
                    81 15
                    82 15
                    83 24
                    84 24
                    85 24
                    86 24
                    87 24
                    88 24
                    89 24
                    90 24
                    91 24
                    92 24
                    93 26
                    94 26
                    95 26
                    96 26
                    97 26
                    98 26
                    99 26
                    100 28
                    101 28
                    102 28
                    103 28
                    104 28
                    105 29
                    106 31
                    107 37
                    108 43
                    109 45
                    110 45
                    111 45
                    112 46
                    113 46
                    114 47
                    115 47
                    116 48
                    117 49
                    118 50
                    119 51
                    120 51
                    121 51
                    122 51
                    123 51
                    124 51
                    125 51
                    126 51
                    127 51
                    128 53
                    129 53
                    130 54
                    131 54
                    132 55
                    133 60
                    134 67
                    135 77
                    136 82
                    137 87
                    138 88
                    139 88
                    140 88
                    141 88
                    142 88
                    143 88
                    144 88
                    145 90
                    146 93
                    147 94
                    148 93
                    149 95
                    150 96
                    151 98
                    152 98
                    153 98
                    154 101
                    155 101
                    156 101
                    157 102
                    158 103
                    159 105
                    160 106
                    161 107
                    162 108
                    163 108
                    164 109
                    165 111
                    166 113
                    167 118
                    168 127
                    169 140
                    170 146
                    171 152
                    172 158
                    173 161
                    174 163
                    175 168
                    176 168
                    177 174
                    178 176
                    179 179
                    180 180
                    181 180
                    182 180
                    183 181
                    184 182
                    185 183
                    186 186
                    187 186
                    188 186
                    189 186
                    190 186
                    191 188
                    192 187
                    193 189
                    194 191
                    195 192
                    196 193
                    197 193
                    198 193
                    199 193
                    200 194
                    201 194
                    202 196
                    203 199
                    204 202
                    205 201
                    206 202
                    207 207
                    208 207
                    209 212
                    210 217
                    211 232
                    212 236
                    213 242
                    214 257
                    215 269
                    216 279
                    217 285
                    218 296
                    219 303
                    220 312
                    221 309
                    222 311
                    223 318
                    224 326
                    225 333
                    226 334
                    227 337
                    228 337
                    229 338
                    230 342
                    231 345
                    232 345
                    233 354
                    234 356
                    235 356
                    236 353
                    237 348
                    238 347
                    239 348
                    240 350
                    241 349
                    242 351];
    end %Properties
    
    methods (Static)
        

        
        
        function numE =convertNucsToEcells(nucVec)
        %% ==================
        %   Name:   findNumE.m
        %   Version:    1.0, Nov. 15th 2013
        %   Author:     Allison Wu, modified by Larry Du
        %   Command:    numE=DuFindNumE(totalNucNum)
        %   Description: Use E cell lineage data from EPIC to find the number of E
        %   cells in the designated nuclei stage.
            
        %% ==================
        % In timingE.mat file, 1st column is time (min), second column is
        % nuclei count, third column     
        %Note, I modified this to accept vectors of nuclei counts

            numE = [];
            for i=(1:length(nucVec))
                for j = (1:length(DuNucsTime.findE))
                    % Find the index for row in divMinNucs
                    % (e.g. 185 nuclei is not listed in table, find nearest
                    % neighbor 184 nuclei)
                    [~,smallestEntryInd]=min(abs(DuNucsTime.findE(:,2)- nucVec(i)));
                    % Find the number of E nuclei at this time point
                    numE(i,1)=DuNucsTime.findE(smallestEntryInd,3);  

                end
            end
        end
        


        function timeMin = convertNucsToTime(nucVec)
            timeMin = [];
            for i=(1:length(nucVec))
                for j = (1:length(DuNucsTime.divTimes))
                    % Find the index for row in divMinNucs
                    % (e.g. 185 nuclei is not listed in table, find nearest
                    % neighbor 184 nuclei)
                    [~,bestInd]=min(abs(DuNucsTime.divTimes(:,2)- nucVec(i)));
                    % Find the time in min corresponding to this number of nuclei
                    timeMin(i,1)=DuNucsTime.divTimes(bestInd,1);  
                    
                end
            end
            
        end %convertNucsToTime
        
    end %Methods
    
end %classdef
