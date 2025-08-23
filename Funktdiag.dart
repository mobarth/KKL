#include <Wire.h>

const int pwmPin = 10;
const int pwmTestValue = 128;
int lastStates[24];
bool pwmAktiv = true;

void setup() {
  Serial.begin(9600);
  while (!Serial);
  Serial.println("🔍 Diagnose gestartet...");

  // PWM-Test
  pinMode(pwmPin, OUTPUT);
  analogWriteFrequency(pwmPin, 25000);
  analogWrite(pwmPin, pwmTestValue);
  Serial.println("⚡ PWM aktiviert auf Pin 10 mit 25kHz, 50% Duty");

  // I2C-Scan
  Wire.begin();
  delay(1000);
  Serial.println("🔎 Starte I2C-Scan...");
  for (byte address = 1; address < 127; address++) {
    Wire.beginTransmission(address);
    if (Wire.endTransmission() == 0) {
      Serial.print("✅ I2C-Gerät gefunden bei Adresse 0x");
      Serial.println(address, HEX);
    }
  }

  for (int pin = 0; pin <= 23; pin++) {
    pinMode(pin, INPUT_PULLUP);
    lastStates[pin] = digitalRead(pin);
  }

  Serial.println("📡 Starte Live-Pin-Überwachung. Sende 'x' zum Stoppen des Lüfters.");
}

void loop() {
  // Tasteneingabe prüfen
  if (Serial.available()) {
    char c = Serial.read();
    if (c == 'x') {
      analogWrite(pwmPin, 0);
      pwmAktiv = false;
      Serial.println("🛑 PWM auf Pin 10 deaktiviert.");
    }
  }

  // Pin-Status prüfen
  for (int pin = 0; pin <= 23; pin++) {
    int state = digitalRead(pin);
    if (state != lastStates[pin]) {
      Serial.print("📍 Pin ");
      Serial.print(pin);
      Serial.print(" wurde ");
      Serial.println(state == HIGH ? "HIGH" : "LOW");
      lastStates[pin] = state;
    }
  }

  delay(200);
}
