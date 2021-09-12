void setup() {
  // put your setup code here, to run once:
  analogWriteResolution (12);
  analogWrite(DAC0,0);
  dacc_set_channel_selection(DACC_INTERFACE, 0);
}

void loop() {
  // put your main code here, to run repeatedly:
  dacc_write_conversion_data(DACC_INTERFACE, 0);
  dacc_write_conversion_data(DACC_INTERFACE, 4095);
}
