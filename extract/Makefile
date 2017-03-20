include ../make.inc
include make_extract.inc

TARGETS=extract_exact_loop_HighMomForm_h5\
	extract_exact_loop_h5\
	extract_stoch_loop_h5\
	extract_stoch_loop_TSM_h5\
	extract_threep_h5\
	extract_threep_posSpace_FastFT_h5\
	extract_twop_baryons_h5\
	extract_twop_baryons_posSpace_FastFT_h5\
	extract_twop_mesons_h5\
	extract_twop_mesons_posSpace_FastFT_h5\
	extract_twop_baryons_HighMomForm_h5\
	extract_twop_mesons_HighMomForm_h5\
	extract_threep_HighMomForm_h5

all: $(TARGETS)

extract_exact_loop_HighMomForm_h5: extract_exact_loop_HighMomForm_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

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

extract_twop_mesons_HighMomForm_h5: extract_twop_mesons_HighMomForm_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

extract_threep_HighMomForm_h5: extract_threep_HighMomForm_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	$(RM) $(TARGETS) ${addsuffix .o, $(TARGETS)}
