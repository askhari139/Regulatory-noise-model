"""
discord_notifier.jl

Sends Discord webhook messages. Reads the webhook URL from the environment
variable DISCORD_WEBHOOK.

All public functions are non-throwing: failures are printed to stderr but
never propagate, so a Discord outage can't kill the pipeline.
"""

using JSON
using Dates

# ============================================================
# CORE SEND
# ============================================================

"""
    send_discord_message(json_body::String)

POST a pre-serialized JSON body to the Discord webhook.
Silently logs errors rather than throwing.
"""
function send_discord_message(json_body::String)
    webhook_url = get(ENV, "DISCORD_WEBHOOK", "")

    if isempty(webhook_url)
        @warn "DISCORD_WEBHOOK not set — skipping notification"
        return
    end

    try
        cmd = `curl -s -o /dev/null -w "%{http_code}" \
                    -H "Content-Type: application/json" \
                    -d $json_body \
                    $webhook_url`
        http_code = strip(read(cmd, String))

        if http_code ∉ ("200", "204")
            @warn "Discord webhook returned HTTP $http_code"
        end
    catch e
        @warn "Discord notification failed: $e"
    end
end

"""
    send_discord_embed(title, description; color, fields)

Build and send a Discord embed message.

# Arguments
- `title`:       Bold heading shown in the embed
- `description`: Body text (supports Discord markdown)
- `color`:       Left-bar colour as an integer (default: blurple 0x5865F2)
- `fields`:      Optional vector of `Dict("name"=>..., "value"=>..., "inline"=>...)` dicts
"""
function send_discord_embed(title::String, description::String;
                             color::Integer=0x5865F2,
                             fields::Vector{Dict{String,Any}}=Dict{String,Any}[])
    payload = Dict(
        "embeds" => [Dict(
            "title"       => title,
            "description" => description,
            "color"       => color,
            "fields"      => fields,
            "timestamp"   => Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"
        )]
    )
    send_discord_message(JSON.json(payload))
end

# convenience alias used inside monitor_and_collect.jl
const discord_embed = send_discord_embed

# ============================================================
# NAMED NOTIFICATION FUNCTIONS  (called from 0_run_all.jl)
# ============================================================

"""Notify that the full pipeline is starting."""
function notify_pipeline_start(networks::Vector{String}, noise_modes::Vector{String})
    combos = ["$nm/$net" for net in networks for nm in noise_modes]
    send_discord_embed(
        "🧬 Pipeline Starting",
        "Running **$(length(combos))** combination(s):\n" *
        join("• " .* combos, "\n"),
        color=0x57F287
    )
end

"""
Notify that the collection step is starting.
Called from monitor_and_collect.jl via notify_collection_starting.
"""
function notify_collection_starting(work_dirs, total_expected, elapsed_secs)
    send_discord_embed(
        "📦 Collecting Results",
        "All **$total_expected** simulations finished in **$(fmt_duration(elapsed_secs))**.\n" *
        "Starting collection + plotting for **$(length(work_dirs))** directory/directories...",
        color=0xEB459E
    )
end

"""Notify that collection is done."""
function notify_collection_complete(n_successful::Int)
    send_discord_embed(
        "✅ Collection Complete",
        "Successfully collected results for **$n_successful** directory/directories.",
        color=0x57F287
    )
end

# ============================================================
# HELPERS (also exported for monitor_and_collect.jl)
# ============================================================

"""Format a duration in seconds to a human-readable string."""
function fmt_duration(secs::Real)
    secs < 60   && return "$(round(Int, secs))s"
    secs < 3600 && return "$(round(secs / 60,   digits=1)) min"
    return "$(round(secs / 3600, digits=2)) h"
end