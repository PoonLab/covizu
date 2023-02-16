module.exports = {
    apps: [{
        name: "app",
        script: "server.js",
        wait_ready: true,
        listen_timeout: 15000,
        instances : 2,
        exec_mode : "cluster"
    }]
};