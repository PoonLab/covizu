module.exports = {
    apps: [{
        name: "app",
        script: "server.js",
        wait_ready: true,
        env_production: {
          NODE_ENV: "PROD"
        },
        listen_timeout: 60000,
        instances : 2,
        exec_mode : "cluster"
    }]
};
