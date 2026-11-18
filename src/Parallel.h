#pragma once

#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace TracEon::detail {

// Join on launch failure; propagate worker failures only after all workers exit.
template <typename Worker>
void parallelFor(size_t begin, size_t end, Worker worker) {
    std::vector<std::exception_ptr> errors(end - begin);
    {
        std::vector<std::jthread> threads;
        threads.reserve(end - begin);
        for (size_t i = begin; i < end; ++i) {
            threads.emplace_back([&, i] {
                try { worker(i); }
                catch (...) { errors[i - begin] = std::current_exception(); }
            });
        }
    }
    for (const auto& error : errors) {
        if (error) std::rethrow_exception(error);
    }
}

} // namespace TracEon::detail
