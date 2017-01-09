include ../make.inc


TARGETS=extract_exact_loop_h5\
	extract_stoch_loop_h5\
	extract_stoch_loop_TSM_h5\
	extract_threep_h5\
	extract_threep_posSpace_FastFT_h5\
	extract_twop_baryons_h5\
	extract_twop_baryons_posSpace_FastFT_h5\
	extract_twop_mesons_h5\
	extract_twop_mesons_posSpace_FastFT_h5\
	extract_twop_baryons_HighMomForm_h5

INCLUDE := -I/opt/cray/pe/fftw/3.3.4.10/haswell/include $(INCLUDE)
LDFLAGS := -L/opt/cray/pe/fftw/3.3.4.10/haswell/lib $(LDFLAGS)
LIBS := -lfftw3f -lfftw3 $(LIBS)


all: $(TARGETS)

extract_exact_loop_h5: extract_exact_loop_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_stoch_loop_h5: extract_stoch_loop_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_stoch_loop_TSM_h5: extract_stoch_loop_TSM_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_threep_h5: extract_threep_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_threep_posSpace_FastFT_h5: extract_threep_posSpace_FastFT_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_twop_baryons_h5: extract_twop_baryons_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_twop_baryons_posSpace_FastFT_h5: extract_twop_baryons_posSpace_FastFT_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_twop_mesons_h5: extract_twop_mesons_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_twop_mesons_posSpace_FastFT_h5: extract_twop_mesons_posSpace_FastFT_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_twop_baryons_HighMomForm_h5: extract_twop_baryons_HighMomForm_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	$(RM) $(TARGETS) ${addsuffix .o, $(TARGETS)}
