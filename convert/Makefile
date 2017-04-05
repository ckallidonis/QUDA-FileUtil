include ../make.inc


TARGETS=convert_exact_loop_h5\
	convert_exact_loop_ultra_local_h5\
	convert_stoch_loop_h5\
	convert_stoch_loop_ultra_local_h5\
	convert_threep_h5\
	convert_threep_multi_source_h5\
	convert_twop_baryons_multi_source_h5\
	convert_twop_mesons_multi_source_h5


all: $(TARGETS)


convert_exact_loop_h5: convert_exact_loop_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_exact_loop_ultra_local_h5: convert_exact_loop_ultra_local_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_stoch_loop_h5: convert_stoch_loop_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_stoch_loop_ultra_local_h5: convert_stoch_loop_ultra_local_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_threep_h5: convert_threep_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_threep_multi_source_h5: convert_threep_multi_source_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_twop_baryons_multi_source_h5: convert_twop_baryons_multi_source_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

convert_twop_mesons_multi_source_h5: convert_twop_mesons_multi_source_h5.o
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)


%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	$(RM) $(TARGETS) ${addsuffix .o, $(TARGETS)}
