// Socle Vulkan compute (phase 1 du plan « Bird entier en C++, puis Vulkan »).
//
// CONTRAT : l'absence de GPU n'est JAMAIS une erreur. Sans APEX_WITH_VULKAN a la
// compilation, sans vulkan-1.dll (charge a l'execution par volk, aucun lien
// statique), sans GPU compatible : `status()` rend ok = false et la RAISON, et
// chaque appelant retombe sur son chemin CPU (meme discipline que
// MeshDistance::resolve_backend).
//
// Cible actuelle : un GPU discret (NVIDIA prefere), une file compute, buffers
// VMA, dispatchs decoupes pour rester sous le seuil de coupure de Windows (TDR,
// ~2 s) si le GPU pilote aussi l'ecran.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace apex::gpu {

struct GpuStatus {
    bool ok = false;
    std::string reason;        // vide si ok
    std::string device;        // nom du GPU retenu
    bool fp64 = false;         // shaderFloat64 disponible ET active
    double timestamp_period_ns = 0.0;
};

// Initialisation paresseuse (une seule fois, sous verrou). Toujours sans exception.
GpuStatus status();

// Vrai si le module a ete compile avec APEX_WITH_VULKAN.
bool compiled_with_vulkan();

// --- Noyaux de base (fumee et micro-bancs) ----------------------------------
// y <- a*x + y, en float32 ou en float64 (selon `fp64`), sans contraction en FMA
// (decoration NoContraction) : resultat identique au bit au calcul CPU.
// Rend le temps GPU mesure par timestamps (secondes), ou -1 si indisponible.
// Leve std::runtime_error si le GPU n'est pas disponible (l'appelant teste status()).
double saxpy_f32(float a, const float* x, float* y, std::int64_t n);
double saxpy_f64(double a, const double* x, double* y, std::int64_t n);

}  // namespace apex::gpu
