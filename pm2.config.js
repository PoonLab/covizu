module.exports = {
    apps: [{
        name: "epicov_app",
        script: "server.js",
        wait_ready: true,
        env_production: {
          NODE_ENV: "PROD"
        },
        listen_timeout: 60000,
        instances : 2,
        exec_mode : "cluster",
        max_memory_restart: "8000M",
        node_args: [
          "--max_old_space_size=8192"
        ]
    }]
};
