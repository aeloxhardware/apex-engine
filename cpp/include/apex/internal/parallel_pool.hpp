// Pool de fils PERSISTANT pour des boucles paralleles tres frequentes (produits
// matrice-vecteur du solveur LSQ de Bird : des dizaines de milliers d'appels de
// l'ordre de la milliseconde, ou creer des fils a chaque appel couterait plus que
// le gain). Structure interne, sans contrat d'equivalence propre : chaque appelant
// garantit que le resultat ne depend pas du decoupage (iterations independantes,
// aucune reduction partagee).
//
// Jamais detruit (fils bloques en attente a la fin du processus) : joindre des
// fils pendant le dechargement d'une DLL Windows peut bloquer (verrou du chargeur).
#pragma once

#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace apex::internal {

class ParallelPool {
public:
    explicit ParallelPool(unsigned threads) : n_(threads < 1 ? 1 : threads) {
        for (unsigned id = 1; id < n_; ++id) workers_.emplace_back([this, id] { loop(id); });
    }
    ParallelPool(const ParallelPool&) = delete;
    ParallelPool& operator=(const ParallelPool&) = delete;

    unsigned size() const { return n_; }

    // fn(begin, end) sur [0, count) decoupe en `size()` blocs contigus ; l'appelant
    // traite le bloc 0. En dessous de `min_parallel`, appel direct (aucune
    // synchronisation). Un seul run() a la fois (verrou d'appel).
    void run(std::int64_t count, const std::function<void(std::int64_t, std::int64_t)>& fn,
             std::int64_t min_parallel = 4096) {
        if (n_ <= 1 || count < min_parallel) {
            fn(0, count);
            return;
        }
        std::lock_guard<std::mutex> call(call_mutex_);
        {
            std::lock_guard<std::mutex> lock(mutex_);
            fn_ = &fn;
            count_ = count;
            remaining_ = n_ - 1;
            ++generation_;
        }
        wake_.notify_all();
        chunk(0, fn);
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [this] { return remaining_ == 0; });
        fn_ = nullptr;
    }

private:
    void chunk(unsigned id, const std::function<void(std::int64_t, std::int64_t)>& fn) const {
        const std::int64_t size = (count_ + n_ - 1) / n_;
        const std::int64_t b = static_cast<std::int64_t>(id) * size;
        const std::int64_t e = b + size < count_ ? b + size : count_;
        if (b < e) fn(b, e);
    }

    void loop(unsigned id) {
        std::uint64_t seen = 0;
        for (;;) {
            const std::function<void(std::int64_t, std::int64_t)>* fn = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait(lock, [&] { return generation_ != seen; });
                seen = generation_;
                fn = fn_;
            }
            chunk(id, *fn);
            {
                std::lock_guard<std::mutex> lock(mutex_);
                if (--remaining_ == 0) done_.notify_one();
            }
        }
    }

    unsigned n_;
    std::vector<std::thread> workers_;
    std::mutex call_mutex_, mutex_;
    std::condition_variable wake_, done_;
    const std::function<void(std::int64_t, std::int64_t)>* fn_ = nullptr;
    std::int64_t count_ = 0;
    unsigned remaining_ = 0;
    std::uint64_t generation_ = 0;
};

}  // namespace apex::internal
