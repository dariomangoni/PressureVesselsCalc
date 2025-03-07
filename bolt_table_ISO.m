% nominal bolt diameter
bolt_table.dB0 = [1, 1.2, 1.4, 1.6, 1.7, 1.8, 2, 2.2, 2.3, 2.5, 2.6, 3, 3.5, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 27, 30, 33, 36, 39, 42, 45, 48, 52, 56, 60, 64, 68, 72, 80, 90, 100];

% bolt thread pitch
bolt_table.P = [0.25, 0.25, 0.3, 0.35, 0.35, 0.35, 0.4, 0.45, 0.4, 0.45, 0.45, 0.5, 0.6, 0.7, 0.8, 1, 1, 1.25, 1.25, 1.5, 1.5, 1.75, 2, 2, 2.5, 2.5, 2.5, 3, 3, 3.5, 3.5, 4, 4, 4.5, 4.5, 5, 5, 5.5, 5.5, 6, 6, 6, 6, 6, 6];

% bolt holes (tap drill size)
bolt_table.dh = bolt_table.dB0 - bolt_table.P;

bolt_table.deltaB = 0.9382*bolt_table.P;

% effective bolt diameter (tensile stress diameter)
bolt_table.dBe = bolt_table.dB0-bolt_table.deltaB;

% shank diameter (cautelative
bolt_table.dBs = bolt_table.dBe;

% cross-section
bolt_table.AB = pi/4*bolt_table.dBe.^2;