from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
from flask_migrate import Migrate
from flask_cors import CORS
from os import path

# Global variables
db = SQLAlchemy()
migrate = Migrate()
DB_NAME = "database.db"

def create_app():
    # Initialize the Flask application
    app = Flask(__name__)
    app.config['SECRET_KEY'] = 'mysecretkeyisolomon'
    app.config['SQLALCHEMY_DATABASE_URI'] = f'sqlite:///{DB_NAME}'
    
    # Initialize extensions
    db.init_app(app)
    migrate.init_app(app, db)
    CORS(app, resources={r"/*": {"origins": "*"}})  # Adjust origins as needed

    # Register Blueprints
    from .routes import routes
    from .auth import auth
    app.register_blueprint(routes, url_prefix='/')
    app.register_blueprint(auth, url_prefix='/')

    # Database setup
    with app.app_context():
        create_database()

    # Login Manager setup
    setup_login_manager(app)

    return app

def create_database():
    """Creates the database if it doesn't exist."""
    if not path.exists('website/' + DB_NAME):
        db.create_all()
        print('Created Database!')

def setup_login_manager(app):
    """Configures the LoginManager for the app."""
    login_manager = LoginManager()
    login_manager.login_view = 'auth.login'
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        from .models import User
        return User.query.get(int(user_id))